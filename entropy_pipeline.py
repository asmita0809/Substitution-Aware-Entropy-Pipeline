#!/usr/bin/env python3
"""
Substitution-Aware Entropy Pipeline  (BLOSUM62)
================================================
Stage 1 — Build entropy profile from MAFFT-aligned natural sequences:
    • Per-column BLOSUM62-weighted substitution probability distribution
    • Per-column Shannon entropy (bits) and normalised entropy [0,1]

Stage 2 — Score designed sequences from a combined MAFFT --add file:
    • Automatically splits natural vs designed by --design-prefix
    • Metric A : Column entropy lookup   — norm_H at that alignment column
    • Metric B : Distribution deviation  — -log2 P(designed AA | BLOSUM natural dist)
    • Combined : mean(A, B)

Typical workflow
----------------
# 1. Align natural sequences
mafft --auto natural.fasta > natural_aligned.fasta

# 2. Add designed sequences into the existing alignment
mafft --add designed.fasta --reorder natural_aligned.fasta > combined_aligned.fasta

# 3. Build natural entropy profile
python entropy_pipeline.py natural \\
    --fasta natural_aligned.fasta \\
    --out natural --plot

# 4. Score designed sequences from the combined file
python entropy_pipeline.py designed \\
    --fasta combined_aligned.fasta \\
    --profile natural_profile.json \\
    --design-prefix "design_" \\
    --out results --plot

# OR run both stages at once
python entropy_pipeline.py both \\
    --natural natural_aligned.fasta \\
    --combined combined_aligned.fasta \\
    --design-prefix "design_" \\
    --out results --plot

Requirements
------------
pip install numpy matplotlib
"""

import argparse
import csv
import json
import math
import sys
from collections import Counter

# ── BLOSUM62 ───────────────────────────────────────────────────────────────
AA_ORDER = list("ARNDCQEGHILKMFPSTWYV")
AA_SET   = set(AA_ORDER)

BLOSUM62 = {
    "A": {"A":4,"R":-1,"N":-2,"D":-2,"C":0,"Q":-1,"E":-1,"G":0,"H":-2,"I":-1,"L":-1,"K":-1,"M":-1,"F":-2,"P":-1,"S":1,"T":0,"W":-3,"Y":-2,"V":0},
    "R": {"A":-1,"R":5,"N":0,"D":-2,"C":-3,"Q":1,"E":0,"G":-2,"H":0,"I":-3,"L":-2,"K":2,"M":-1,"F":-3,"P":-2,"S":-1,"T":-1,"W":-3,"Y":-2,"V":-3},
    "N": {"A":-2,"R":0,"N":6,"D":1,"C":-3,"Q":0,"E":0,"G":0,"H":1,"I":-3,"L":-3,"K":0,"M":-2,"F":-3,"P":-2,"S":1,"T":0,"W":-4,"Y":-2,"V":-3},
    "D": {"A":-2,"R":-2,"N":1,"D":6,"C":-3,"Q":0,"E":2,"G":-1,"H":-1,"I":-3,"L":-4,"K":-1,"M":-3,"F":-3,"P":-1,"S":0,"T":-1,"W":-4,"Y":-3,"V":-3},
    "C": {"A":0,"R":-3,"N":-3,"D":-3,"C":9,"Q":-3,"E":-4,"G":-3,"H":-3,"I":-1,"L":-1,"K":-3,"M":-1,"F":-2,"P":-3,"S":-1,"T":-1,"W":-2,"Y":-2,"V":-1},
    "Q": {"A":-1,"R":1,"N":0,"D":0,"C":-3,"Q":5,"E":2,"G":-2,"H":0,"I":-3,"L":-2,"K":1,"M":0,"F":-3,"P":-1,"S":0,"T":-1,"W":-2,"Y":-1,"V":-2},
    "E": {"A":-1,"R":0,"N":0,"D":2,"C":-4,"Q":2,"E":5,"G":-2,"H":0,"I":-3,"L":-3,"K":1,"M":-2,"F":-3,"P":-1,"S":0,"T":-1,"W":-3,"Y":-2,"V":-2},
    "G": {"A":0,"R":-2,"N":0,"D":-1,"C":-3,"Q":-2,"E":-2,"G":6,"H":-2,"I":-4,"L":-4,"K":-2,"M":-3,"F":-3,"P":-2,"S":0,"T":-2,"W":-2,"Y":-3,"V":-3},
    "H": {"A":-2,"R":0,"N":1,"D":-1,"C":-3,"Q":0,"E":0,"G":-2,"H":8,"I":-3,"L":-3,"K":-1,"M":-2,"F":-1,"P":-2,"S":-1,"T":-2,"W":-2,"Y":2,"V":-3},
    "I": {"A":-1,"R":-3,"N":-3,"D":-3,"C":-1,"Q":-3,"E":-3,"G":-4,"H":-3,"I":4,"L":2,"K":-3,"M":1,"F":0,"P":-3,"S":-2,"T":-1,"W":-3,"Y":-1,"V":3},
    "L": {"A":-1,"R":-2,"N":-3,"D":-4,"C":-1,"Q":-2,"E":-3,"G":-4,"H":-3,"I":2,"L":4,"K":-2,"M":2,"F":0,"P":-3,"S":-2,"T":-1,"W":-2,"Y":-1,"V":1},
    "K": {"A":-1,"R":2,"N":0,"D":-1,"C":-3,"Q":1,"E":1,"G":-2,"H":-1,"I":-3,"L":-2,"K":5,"M":-1,"F":-3,"P":-1,"S":0,"T":-1,"W":-3,"Y":-2,"V":-2},
    "M": {"A":-1,"R":-1,"N":-2,"D":-3,"C":-1,"Q":0,"E":-2,"G":-3,"H":-2,"I":1,"L":2,"K":-1,"M":5,"F":0,"P":-2,"S":-1,"T":-1,"W":-1,"Y":-1,"V":1},
    "F": {"A":-2,"R":-3,"N":-3,"D":-3,"C":-2,"Q":-3,"E":-3,"G":-3,"H":-1,"I":0,"L":0,"K":-3,"M":0,"F":6,"P":-4,"S":-2,"T":-2,"W":1,"Y":3,"V":-1},
    "P": {"A":-1,"R":-2,"N":-2,"D":-1,"C":-3,"Q":-1,"E":-1,"G":-2,"H":-2,"I":-3,"L":-3,"K":-1,"M":-2,"F":-4,"P":7,"S":-1,"T":-1,"W":-4,"Y":-3,"V":-2},
    "S": {"A":1,"R":-1,"N":1,"D":0,"C":-1,"Q":0,"E":0,"G":0,"H":-1,"I":-2,"L":-2,"K":0,"M":-1,"F":-2,"P":-1,"S":4,"T":1,"W":-3,"Y":-2,"V":-2},
    "T": {"A":0,"R":-1,"N":0,"D":-1,"C":-1,"Q":-1,"E":-1,"G":-2,"H":-2,"I":-1,"L":-1,"K":-1,"M":-1,"F":-2,"P":-1,"S":1,"T":5,"W":-2,"Y":-2,"V":0},
    "W": {"A":-3,"R":-3,"N":-4,"D":-4,"C":-2,"Q":-2,"E":-3,"G":-2,"H":-2,"I":-3,"L":-2,"K":-3,"M":-1,"F":1,"P":-4,"S":-3,"T":-2,"W":11,"Y":2,"V":-3},
    "Y": {"A":-2,"R":-2,"N":-2,"D":-3,"C":-2,"Q":-1,"E":-2,"G":-3,"H":2,"I":-1,"L":-1,"K":-2,"M":-1,"F":3,"P":-3,"S":-2,"T":-2,"W":2,"Y":7,"V":-1},
    "V": {"A":0,"R":-3,"N":-3,"D":-3,"C":-1,"Q":-2,"E":-2,"G":-3,"H":-3,"I":3,"L":1,"K":-2,"M":1,"F":-1,"P":-2,"S":-2,"T":0,"W":-3,"Y":-1,"V":4},
}

MAX_ENTROPY = math.log2(20)   # ~4.3219 bits
GAP_CHARS   = set("-. ")


# ══════════════════════════════════════════════════════════════════════════
# I/O helpers
# ══════════════════════════════════════════════════════════════════════════

def read_fasta(path: str) -> list[tuple[str, str]]:
    """Return list of (header, sequence) from a FASTA file."""
    records, header, seq_parts = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header, seq_parts = line[1:], []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records


def split_combined(
    records: list[tuple[str, str]],
    design_prefix: str,
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """
    Split a combined FASTA (natural + designed from mafft --add) into
    (natural_records, designed_records) by matching the header prefix.

    Matching is case-insensitive so 'Design_' matches 'design_001'.
    The prefix is compared against the FULL header string (everything
    after the '>' character, including description fields).
    """
    prefix  = design_prefix.lower()
    natural, designed = [], []
    for header, seq in records:
        if header.lower().startswith(prefix):
            designed.append((header, seq))
        else:
            natural.append((header, seq))
    return natural, designed


# ══════════════════════════════════════════════════════════════════════════
# Stage 1 — Natural entropy profile
# ══════════════════════════════════════════════════════════════════════════

def blosum62_softmax(aa: str) -> list[float]:
    """BLOSUM62 row → substitution probability via softmax."""
    row    = BLOSUM62[aa]
    scores = [row[b] for b in AA_ORDER]
    max_s  = max(scores)
    exps   = [math.exp(s - max_s) for s in scores]
    total  = sum(exps)
    return [e / total for e in exps]


def column_blosum62_distribution(column_aas: list[str]) -> list[float]:
    """
    Frequency-weighted mixture of BLOSUM62 softmax vectors for one column.
    Produces a substitution-smoothed distribution over all 20 AAs.
    """
    valid = [aa for aa in column_aas if aa in AA_SET]
    if not valid:
        return [1.0 / 20] * 20          # all-gap column → uniform
    counts = Counter(valid)
    total  = sum(counts.values())
    mixed  = [0.0] * 20
    for aa, cnt in counts.items():
        w     = cnt / total
        probs = blosum62_softmax(aa)
        mixed = [m + w * p for m, p in zip(mixed, probs)]
    return mixed


def shannon_entropy(probs: list[float]) -> float:
    return -sum(p * math.log2(p) for p in probs if p > 1e-12)


def compute_natural_profile(records: list[tuple[str, str]]) -> list[dict]:
    """
    Compute per-column entropy profile from aligned natural sequences.
    """
    if not records:
        raise ValueError("No natural sequences found.")

    aln_len = len(records[0][1])
    seqs    = [s for _, s in records]
    profile = []

    for col in range(aln_len):
        col_chars   = [s[col] for s in seqs if col < len(s)]
        gap_count   = sum(1 for c in col_chars if c in GAP_CHARS)
        gap_frac    = gap_count / len(col_chars) if col_chars else 1.0
        col_aas     = [c for c in col_chars if c not in GAP_CHARS]

        freq_counts = Counter(c for c in col_aas if c in AA_SET)
        n_valid     = sum(freq_counts.values())
        freq_dist   = (
            {aa: freq_counts.get(aa, 0) / n_valid for aa in AA_ORDER}
            if n_valid else {aa: 0.05 for aa in AA_ORDER}
        )

        mixed_dist = column_blosum62_distribution(col_aas)
        H          = shannon_entropy(mixed_dist)
        norm_H     = H / MAX_ENTROPY

        profile.append({
            "col":         col,
            "position":    col + 1,
            "gap_frac":    round(gap_frac,  6),
            "n_seqs":      len(col_aas),
            "freq_dist":   freq_dist,
            "blosum_dist": {aa: round(p, 8) for aa, p in zip(AA_ORDER, mixed_dist)},
            "H_bits":      round(H,      8),
            "norm_H":      round(norm_H, 8),
        })

    return profile


def save_natural_profile(profile: list[dict], prefix: str) -> str:
    json_path = f"{prefix}_profile.json"
    with open(json_path, "w") as f:
        json.dump(profile, f, indent=2)
    print(f"  ✓ Natural entropy profile (JSON) → {json_path}")

    csv_path = f"{prefix}_profile.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            ["position", "gap_frac", "n_seqs", "H_bits", "norm_H"] +
            [f"freq_{aa}"   for aa in AA_ORDER] +
            [f"blosum_{aa}" for aa in AA_ORDER]
        )
        for c in profile:
            w.writerow(
                [c["position"], c["gap_frac"], c["n_seqs"], c["H_bits"], c["norm_H"]] +
                [round(c["freq_dist"].get(aa,   0), 6) for aa in AA_ORDER] +
                [round(c["blosum_dist"].get(aa, 0), 8) for aa in AA_ORDER]
            )
    print(f"  ✓ Natural entropy profile (CSV)  → {csv_path}")
    return json_path


# ══════════════════════════════════════════════════════════════════════════
# Stage 2 — Designed sequence penalty scoring
# ══════════════════════════════════════════════════════════════════════════

def score_designed_sequence(seq: str, profile: list[dict]) -> list[dict]:
    """
    Score every non-gap residue of a designed sequence.

    Metric A — entropy lookup:
        penalty_A = norm_H of that alignment column
        Low  → conserved position in nature  (high constraint on design)
        High → variable position in nature   (low constraint)

    Metric B — distribution deviation:
        penalty_B = -log2( P(designed_AA | BLOSUM62-smoothed natural dist) )
                    normalised by log2(20) → [0, ~1]
        Low  → designed AA well-supported by natural substitution patterns
        High → designed AA is unexpected / unusual

    Combined = mean(A, B)
    """
    scores = []
    for col, prof in enumerate(profile):
        if col >= len(seq):
            break
        aa = seq[col]
        if aa in GAP_CHARS:
            continue

        aa     = aa.upper()
        norm_H = prof["norm_H"]

        p_aa          = max(prof["blosum_dist"].get(aa, 1e-8), 1e-8)
        neg_log_p     = -math.log2(p_aa)
        norm_neg_logp = neg_log_p / MAX_ENTROPY
        combined      = (norm_H + norm_neg_logp) / 2.0

        scores.append({
            "col":                        col,
            "position":                   col + 1,
            "aa":                         aa,
            "gap_frac_natural":           prof["gap_frac"],
            "H_bits":                     prof["H_bits"],
            "penalty_A_entropy_lookup":   round(norm_H,          6),
            "penalty_B_distribution_dev": round(norm_neg_logp,   6),
            "combined_penalty":           round(combined,         6),
            "p_aa_in_natural":            round(p_aa,             8),
        })
    return scores


def aggregate_scores(scores: list[dict]) -> dict:
    if not scores:
        return {}
    keys = [
        "penalty_A_entropy_lookup",
        "penalty_B_distribution_dev",
        "combined_penalty",
    ]
    return {
        "n_positions": len(scores),
        **{f"mean_{k}": round(sum(s[k] for s in scores) / len(scores), 6) for k in keys},
        **{f"sum_{k}":  round(sum(s[k] for s in scores),               6) for k in keys},
        **{f"max_{k}":  round(max(s[k] for s in scores),               6) for k in keys},
        **{f"min_{k}":  round(min(s[k] for s in scores),               6) for k in keys},
    }


def score_sequence_group(
    records: list[tuple[str, str]],
    profile: list[dict],
    prefix:  str,
    group:   str,
) -> list[dict]:
    group_slug = group.lower().strip().replace(" ", "_")
    all_results = []

    for header, seq in records:
        scores = score_designed_sequence(seq, profile)
        agg    = aggregate_scores(scores)
        all_results.append({
            "header":       header,
            "sequence":     seq.replace("-", "").replace(".", ""),
            "aggregate":    agg,
            "per_position": scores,
        })

    # per-sequence JSON
    json_path = f"{prefix}_{group_slug}_scores.json"
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"  ✓ {group.title()} scores (JSON)            → {json_path}")

    # per-sequence summary CSV
    if all_results:
        sum_csv = f"{prefix}_{group_slug}_summary.csv"
        summary_rows = [{"sequence": r["header"], **r["aggregate"]} for r in all_results]
        with open(sum_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
            w.writeheader()
            w.writerows(summary_rows)
        print(f"  ✓ {group.title()} summary CSV              → {sum_csv}")

    # long-format per-position CSV
    pos_csv = f"{prefix}_{group_slug}_per_position.csv"
    with open(pos_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "sequence", "col", "position", "aa",
            "H_bits", "gap_frac_natural",
            "penalty_A_entropy_lookup",
            "penalty_B_distribution_dev",
            "combined_penalty",
            "p_aa_in_natural",
        ])
        for res in all_results:
            for s in res["per_position"]:
                w.writerow([
                    res["header"],
                    s["col"], s["position"], s["aa"],
                    s["H_bits"], s["gap_frac_natural"],
                    s["penalty_A_entropy_lookup"],
                    s["penalty_B_distribution_dev"],
                    s["combined_penalty"],
                    s["p_aa_in_natural"],
                ])
    print(f"  ✓ {group.title()} per-position CSV         → {pos_csv}")

    return all_results


def score_all_designed(
    records: list[tuple[str, str]],
    profile: list[dict],
    prefix:  str,
) -> list[dict]:
    return score_sequence_group(records, profile, prefix, "designed")


def score_all_natural(
    records: list[tuple[str, str]],
    profile: list[dict],
    prefix:  str,
) -> list[dict]:
    return score_sequence_group(records, profile, prefix, "natural")


# ══════════════════════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════════════════════

def _style_axes(axes):
    for ax in axes:
        ax.set_facecolor("#111827")
        for sp in ax.spines.values():
            sp.set_edgecolor("#374151")
        ax.tick_params(colors="#9ca3af", labelsize=8)
        ax.yaxis.label.set_color("#9ca3af")
        ax.xaxis.label.set_color("#9ca3af")
        ax.title.set_color("#e2e8f0")


def _entropy_cmap():
    import matplotlib.colors as mc
    return mc.LinearSegmentedColormap.from_list(
        "ep", [(0,"#14b8a6"),(0.33,"#eab308"),(0.66,"#f97316"),(1,"#dc2626")]
    )


def plot_natural_profile(profile: list[dict], prefix: str):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  ⚠  pip install matplotlib", file=sys.stderr); return

    positions = [c["position"] for c in profile]
    norm_H    = [c["norm_H"]   for c in profile]
    H_bits    = [c["H_bits"]   for c in profile]
    gap_frac  = [c["gap_frac"] for c in profile]
    cmap      = _entropy_cmap()

    fig, axes = plt.subplots(3, 1,
                             figsize=(max(12, len(positions) * 0.12), 10),
                             facecolor="#0b0f1a")
    fig.suptitle("Natural Sequence Entropy Profile  ·  BLOSUM62-Weighted",
                 color="#e2e8f0", fontsize=12, fontweight="bold")
    _style_axes(axes)

    axes[0].bar(positions, norm_H, color=[cmap(v) for v in norm_H], width=1.0)
    axes[0].axhline(sum(norm_H)/len(norm_H), color="#6b7280", lw=1, ls="--",
                    label=f"mean = {sum(norm_H)/len(norm_H):.3f}")
    axes[0].set_ylabel("Normalised Entropy")
    axes[0].set_title("Per-Column Normalised Entropy  (0 = conserved  ·  1 = variable)")
    axes[0].legend(fontsize=8, labelcolor="#9ca3af", facecolor="#1f2937", edgecolor="#374151")

    axes[1].plot(positions, H_bits, color="#14b8a6", lw=1.2)
    axes[1].axhline(MAX_ENTROPY, color="#ef4444", lw=1, ls=":",
                    label=f"max H = {MAX_ENTROPY:.3f} bits")
    axes[1].axhline(sum(H_bits)/len(H_bits), color="#6b7280", lw=1, ls="--",
                    label=f"mean H = {sum(H_bits)/len(H_bits):.3f} bits")
    axes[1].set_ylabel("Shannon Entropy (bits)")
    axes[1].set_title("Per-Column Shannon Entropy")
    axes[1].legend(fontsize=8, labelcolor="#9ca3af", facecolor="#1f2937", edgecolor="#374151")

    axes[2].fill_between(positions, gap_frac, color="#374151", alpha=0.9)
    axes[2].set_ylabel("Gap Fraction")
    axes[2].set_title("Column Gap Fraction")
    axes[2].set_xlabel("Alignment Column")
    axes[2].set_ylim(0, 1)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = f"{prefix}_natural_profile.png"
    plt.savefig(out, dpi=150, bbox_inches="tight", facecolor="#0b0f1a")
    plt.close()
    print(f"  ✓ Natural profile plot              → {out}")


def plot_designed_scores(
    all_results: list[dict],
    profile:     list[dict],
    prefix:      str,
):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("  ⚠  pip install matplotlib numpy", file=sys.stderr); return

    cmap  = _entropy_cmap()
    names = [r["header"][:30] for r in all_results]
    x     = list(range(len(names)))

    def _get(r, key): return r["aggregate"].get(key, 0)
    meanA    = [_get(r, "mean_penalty_A_entropy_lookup")   for r in all_results]
    meanB    = [_get(r, "mean_penalty_B_distribution_dev") for r in all_results]
    combined = [_get(r, "mean_combined_penalty")           for r in all_results]

    # ── bar chart summary ──
    fig, ax = plt.subplots(figsize=(max(8, len(names)*0.7+2), 5), facecolor="#0b0f1a")
    _style_axes([ax])
    bw = 0.28
    ax.bar([i-bw for i in x], meanA,    width=bw, color="#14b8a6", alpha=0.85,
           label="Penalty A  (entropy lookup)")
    ax.bar([i    for i in x], meanB,    width=bw, color="#f97316", alpha=0.85,
           label="Penalty B  (distribution deviation)")
    ax.bar([i+bw for i in x], combined, width=bw, color="#a855f7", alpha=0.85,
           label="Combined")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean Penalty")
    ax.set_title("Designed Sequence Mean Penalties")
    ax.legend(fontsize=8, labelcolor="#9ca3af", facecolor="#1f2937", edgecolor="#374151")
    plt.tight_layout()
    out1 = f"{prefix}_designed_summary.png"
    plt.savefig(out1, dpi=150, bbox_inches="tight", facecolor="#0b0f1a")
    plt.close()
    print(f"  ✓ Designed summary bar chart        → {out1}")

    # ── heatmap: designed sequences × alignment columns ──
    n_cols = len(profile)
    n_seqs = len(all_results)
    mat    = np.full((n_seqs, n_cols), np.nan)
    for si, res in enumerate(all_results):
        for pos in res["per_position"]:
            mat[si, pos["col"]] = pos["combined_penalty"]

    # natural entropy reference row included below the heatmap
    nat_row = np.array([[c["norm_H"] for c in profile]])

    fig, axes = plt.subplots(
        2, 1,
        figsize=(max(14, n_cols*0.15), max(4, n_seqs*0.45+3)),
        facecolor="#0b0f1a",
        gridspec_kw={"height_ratios": [n_seqs, 1]},
    )
    _style_axes(axes)

    im = axes[0].imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                        interpolation="nearest")
    axes[0].set_yticks(range(n_seqs))
    axes[0].set_yticklabels([r["header"][:30] for r in all_results], fontsize=7)
    axes[0].set_title("Per-Position Combined Penalty  (designed sequences)")
    axes[0].set_xlabel("Alignment Column")
    cbar = fig.colorbar(im, ax=axes[0], fraction=0.015, pad=0.02)
    cbar.ax.tick_params(colors="#9ca3af", labelsize=8)
    cbar.set_label("Combined Penalty", color="#9ca3af", fontsize=8)

    axes[1].imshow(nat_row, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                   interpolation="nearest")
    axes[1].set_yticks([0])
    axes[1].set_yticklabels(["Natural\nEntropy"], fontsize=7)
    axes[1].set_title("Natural Column Entropy  (reference)")
    axes[1].set_xlabel("Alignment Column")

    plt.tight_layout()
    out2 = f"{prefix}_designed_heatmap.png"
    plt.savefig(out2, dpi=150, bbox_inches="tight", facecolor="#0b0f1a")
    plt.close()
    print(f"  ✓ Designed heatmap (+ natural ref)  → {out2}")


# ══════════════════════════════════════════════════════════════════════════
# Terminal print helpers
# ══════════════════════════════════════════════════════════════════════════

def print_natural_summary(profile: list[dict], n_seqs: int):
    norm_H = [c["norm_H"] for c in profile]
    H_bits = [c["H_bits"] for c in profile]
    print(f"\n  {'─'*55}")
    print(f"  Stage 1 — Natural Entropy Profile")
    print(f"  {'─'*55}")
    print(f"  Sequences analysed  : {n_seqs}")
    print(f"  Alignment columns   : {len(profile)}")
    print(f"  Mean norm entropy   : {sum(norm_H)/len(norm_H):.6f}")
    print(f"  Max  norm entropy   : {max(norm_H):.6f}")
    print(f"  Min  norm entropy   : {min(norm_H):.6f}")
    print(f"  Mean H (bits)       : {sum(H_bits)/len(H_bits):.6f}")
    print(f"  {'─'*55}\n")


def print_designed_summary(all_results: list[dict]):
    print(f"\n  {'─'*68}")
    print(f"  Stage 2 — Designed Sequence Penalties")
    print(f"  {'─'*68}")
    print(f"  {'Sequence':<36} {'Pen-A':>8} {'Pen-B':>8} {'Combined':>9}")
    print(f"  {'─'*35} {'─'*7} {'─'*7} {'─'*8}")
    for r in all_results:
        agg  = r["aggregate"]
        name = r["header"][:35]
        pA   = agg.get("mean_penalty_A_entropy_lookup",   0)
        pB   = agg.get("mean_penalty_B_distribution_dev", 0)
        pc   = agg.get("mean_combined_penalty",           0)
        print(f"  {name:<36} {pA:>8.4f} {pB:>8.4f} {pc:>9.4f}")
    print(f"  {'─'*68}")
    print(f"  Penalty A  = normalised column entropy  (natural variability at position)")
    print(f"  Penalty B  = -log2 P(designed AA | BLOSUM62 natural dist)")
    print(f"  Combined   = mean(A, B)")
    print(f"  All values in [0, 1].  Lower = more natural-like.\n")


def print_natural_penalty_summary(all_results: list[dict]):
    print(f"\n  {'─'*68}")
    print(f"  Stage 2B — Natural Sequence Penalties")
    print(f"  {'─'*68}")
    print(f"  {'Sequence':<36} {'Pen-A':>8} {'Pen-B':>8} {'Combined':>9}")
    print(f"  {'─'*35} {'─'*7} {'─'*7} {'─'*8}")
    for r in all_results:
        agg  = r["aggregate"]
        name = r["header"][:35]
        pA   = agg.get("mean_penalty_A_entropy_lookup",   0)
        pB   = agg.get("mean_penalty_B_distribution_dev", 0)
        pc   = agg.get("mean_combined_penalty",           0)
        print(f"  {name:<36} {pA:>8.4f} {pB:>8.4f} {pc:>9.4f}")
    print(f"  {'─'*68}")
    print(f"  Same interpretation as designed scoring; this serves as a natural baseline.\n")


def print_split_report(natural, designed, prefix):
    print(f"  Sequences matching prefix '{prefix}' → scored as designed : {len(designed)}")
    print(f"  Remaining sequences (natural, ignored)                    : {len(natural)}")
    if not designed:
        print(f"\n  ✗ ERROR: No sequences found with prefix '{prefix}'.")
        print(f"    First 8 headers seen in file:")
        for h, _ in (natural + designed)[:8]:
            print(f"      {h}")
        sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════
# CLI commands
# ══════════════════════════════════════════════════════════════════════════

def cmd_natural(args):
    print(f"\n[Stage 1] Natural sequences: {args.fasta}")
    records = read_fasta(args.fasta)
    print(f"  Loaded {len(records)} sequences  |  alignment length {len(records[0][1])}")
    profile = compute_natural_profile(records)
    print_natural_summary(profile, len(records))
    save_natural_profile(profile, args.out)
    if args.plot:
        plot_natural_profile(profile, args.out)


def cmd_designed(args):
    print(f"\n[Stage 2] Natural profile  : {args.profile}")
    with open(args.profile) as f:
        profile = json.load(f)
    print(f"  Profile columns: {len(profile)}")

    print(f"\n[Stage 2] Combined alignment: {args.fasta}")
    all_records = read_fasta(args.fasta)
    print(f"  Total sequences in file: {len(all_records)}")

    natural, designed = split_combined(all_records, args.design_prefix)
    print_split_report(natural, designed, args.design_prefix)

    all_results = score_all_designed(designed, profile, args.out)
    natural_results = score_all_natural(natural, profile, args.out)
    print_designed_summary(all_results)
    print_natural_penalty_summary(natural_results)
    if args.plot:
        plot_designed_scores(all_results, profile, args.out)


def cmd_both(args):
    # Stage 1
    print(f"\n[Stage 1] Natural sequences: {args.natural}")
    nat_records = read_fasta(args.natural)
    print(f"  Loaded {len(nat_records)} sequences  |  alignment length {len(nat_records[0][1])}")
    profile = compute_natural_profile(nat_records)
    print_natural_summary(profile, len(nat_records))
    save_natural_profile(profile, args.out)
    if args.plot:
        plot_natural_profile(profile, args.out)

    # Stage 2
    print(f"\n[Stage 2] Combined alignment: {args.combined}")
    all_records = read_fasta(args.combined)
    print(f"  Total sequences in file: {len(all_records)}")

    natural, designed = split_combined(all_records, args.design_prefix)
    print_split_report(natural, designed, args.design_prefix)

    all_results = score_all_designed(designed, profile, args.out)
    natural_results = score_all_natural(natural, profile, args.out)
    print_designed_summary(all_results)
    print_natural_penalty_summary(natural_results)
    if args.plot:
        plot_designed_scores(all_results, profile, args.out)


# ══════════════════════════════════════════════════════════════════════════
# Entry point
# ══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Substitution-Aware Entropy Pipeline  (BLOSUM62)",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # ── natural ──
    p1 = sub.add_parser("natural",
        help="Stage 1: build entropy profile from aligned natural sequences")
    p1.add_argument("--fasta", required=True,
        help="MAFFT-aligned natural sequences (.fasta / .fa)")
    p1.add_argument("--out",   default="natural",
        help="Output file prefix  (default: natural)")
    p1.add_argument("--plot",  action="store_true",
        help="Save profile plots as PNG")

    # ── designed ──
    p2 = sub.add_parser("designed",
        help="Stage 2: score designed sequences from a combined mafft --add file")
    p2.add_argument("--fasta",          required=True,
        help="Combined aligned FASTA produced by:\n"
             "  mafft --add designed.fasta --reorder natural_aligned.fasta > combined.fasta")
    p2.add_argument("--profile",        required=True,
        help="JSON profile produced by Stage 1  (*_profile.json)")
    p2.add_argument("--design-prefix",  required=True,
        help="Header prefix that identifies designed sequences (case-insensitive)\n"
             "Examples:\n"
             "  --design-prefix 'design_'\n"
             "  --design-prefix 'halluc'\n"
             "  --design-prefix 'AF_'")
    p2.add_argument("--out",            default="results",
        help="Output file prefix  (default: results)")
    p2.add_argument("--plot",           action="store_true",
        help="Save score plots as PNG")

    # ── both ──
    p3 = sub.add_parser("both",
        help="Run Stage 1 + Stage 2 in one command")
    p3.add_argument("--natural",        required=True,
        help="MAFFT-aligned natural sequences (.fasta / .fa)")
    p3.add_argument("--combined",       required=True,
        help="Combined FASTA from mafft --add (natural + designed)")
    p3.add_argument("--design-prefix",  required=True,
        help="Header prefix that identifies designed sequences (case-insensitive)")
    p3.add_argument("--out",            default="results",
        help="Output file prefix  (default: results)")
    p3.add_argument("--plot",           action="store_true",
        help="Save all plots as PNG")

    args = parser.parse_args()
    {"natural": cmd_natural, "designed": cmd_designed, "both": cmd_both}[args.cmd](args)
    print("Done.\n")


if __name__ == "__main__":
    main()
