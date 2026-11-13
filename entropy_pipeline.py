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
from pathlib import Path

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


def compute_natural_profile(
    records:      list[tuple[str, str]],
    max_gap_frac: float = 1.0,
) -> list[dict]:
    """
    Compute per-column entropy profile from aligned natural sequences.

    max_gap_frac marks columns as unreliable rather than removing them:
    a column where most natural sequences carry a gap has its
    distribution estimated from a handful of residues, and will happily
    emit the largest penalties in a run from pure sampling noise. Such
    columns get scorable=False and are skipped during scoring.
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
            "scorable":    bool(gap_frac <= max_gap_frac),
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
            ["position", "gap_frac", "n_seqs", "scorable", "H_bits", "norm_H"] +
            [f"freq_{aa}"   for aa in AA_ORDER] +
            [f"blosum_{aa}" for aa in AA_ORDER]
        )
        for c in profile:
            w.writerow(
                [c["position"], c["gap_frac"], c["n_seqs"],
                 int(c.get("scorable", True)), c["H_bits"], c["norm_H"]] +
                [round(c["freq_dist"].get(aa,   0), 6) for aa in AA_ORDER] +
                [round(c["blosum_dist"].get(aa, 0), 8) for aa in AA_ORDER]
            )
    print(f"  ✓ Natural entropy profile (CSV)  → {csv_path}")
    return json_path


# ══════════════════════════════════════════════════════════════════════════
# Column selection — gap filtering and designable-position subsetting
# ══════════════════════════════════════════════════════════════════════════

def parse_positions(
    spec:   str | None,
    n_cols: int,
    tmap:   "TemplateMap | None" = None,
) -> set[int] | None:
    """
    Parse a --positions argument into a set of 0-based column indices.

    Accepted forms:
      None / "all"     → None  (score every column)
      "auto"           → None here; resolved from the designed records,
                         which is the usual choice for a design set
      "80,82,94-99"    → explicit 1-based ALIGNMENT positions
      path/to/file.txt → one 1-based position per line, '#' comments ok
      path/to/x.resfile→ a Rosetta resfile; its designable positions are
                         mapped onto alignment columns via `tmap`, which
                         therefore has to be supplied (--template)

    Bare numbers are alignment column numbers, not residue numbers. A
    resfile is in the template's residue numbering, so it is only usable
    once a TemplateMap pins the two together.
    """
    if spec is None:
        return None
    key = spec.strip().lower()
    if key in ("auto", "all", ""):
        return None

    text = spec
    path = Path(spec)

    if path.exists() and looks_like_resfile(path):
        if tmap is None:
            raise ValueError(
                f"{spec} looks like a Rosetta resfile, which is written in "
                "residue numbering. Pass --template (the template's header in "
                "the alignment) so the residue numbers can be mapped onto "
                "alignment columns."
            )
        rows = parse_resfile(path)
        designable = {
            r for r, (_chain, cmd, _aas) in rows.items()
            if cmd in ("PIKAA", "NOTAA", "ALLAA", "ALLAAwc", "ALLAAxc",
                       "POLAR", "APOLAR")
        }
        cols    = {tmap.res_to_col[r] for r in designable if r in tmap.res_to_col}
        missing = sorted(r for r in designable if r not in tmap.res_to_col)
        print(f"  Resfile {path.name}: {len(rows)} entries, "
              f"{len(designable)} designable → {len(cols)} alignment columns")
        if missing:
            lo, hi = tmap.res_range
            print(f"    !! {len(missing)} residue numbers are outside the "
                  f"template ({lo}-{hi}): {missing[:8]}")
            print( "       Check --template-start against the PDB numbering.")
        if not cols:
            raise ValueError(
                "No resfile position mapped onto the template. The numbering "
                "offset is almost certainly wrong — check --template-start."
            )
        return cols

    if path.exists():
        text = ",".join(
            line.split()[0]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )

    cols: set[int] = set()
    for chunk in text.replace(" ", "").split(","):
        if not chunk:
            continue
        if "-" in chunk.lstrip("-"):
            lo, hi = chunk.split("-", 1)
            cols.update(range(int(lo) - 1, int(hi)))
        else:
            cols.add(int(chunk) - 1)

    bad = sorted(c + 1 for c in cols if not 0 <= c < n_cols)
    if bad:
        raise ValueError(
            f"--positions lists positions outside the alignment (1-{n_cols}): "
            f"{bad[:10]}"
        )
    return cols


def detect_variable_columns(records: list[tuple[str, str]]) -> set[int]:
    """
    Columns at which the supplied sequences differ from one another.

    For a fixed-backbone design set this recovers the designable positions.
    Every other column holds the same residue in every sequence, so it adds
    an identical constant to each per-sequence mean — padding the average
    and shrinking the spread between designs without carrying any signal.
    """
    if not records:
        return set()
    seqs = [s for _, s in records]
    n    = min(len(s) for s in seqs)
    return {col for col in range(n) if len({s[col] for s in seqs}) > 1}


def resolve_scoring_columns(
    profile:   list[dict],
    positions: set[int] | None,
    quiet:     bool = False,
) -> list[int]:
    """
    Intersect the requested subset with the columns the profile considers
    reliable, and say out loud what was dropped.
    """
    n_cols   = len(profile)
    wanted   = set(range(n_cols)) if positions is None else set(positions)
    scorable = {c for c in wanted if profile[c].get("scorable", True)}
    dropped  = sorted(wanted - scorable)

    if not quiet:
        print(f"  Scoring columns     : {len(scorable)} of {n_cols}")
        if positions is not None:
            print(f"    requested subset  : {len(wanted)}")
        if dropped:
            shown = ", ".join(str(profile[c]["position"]) for c in dropped[:12])
            more  = f"  +{len(dropped) - 12} more" if len(dropped) > 12 else ""
            print(f"    dropped as gappy  : {len(dropped)}  [{shown}{more}]")

    if not scorable:
        raise ValueError(
            "No scorable columns remain. Raise --max-gap-frac or widen --positions."
        )
    return sorted(scorable)


# ══════════════════════════════════════════════════════════════════════════
# Stage 2 — Designed sequence penalty scoring
# ══════════════════════════════════════════════════════════════════════════

def score_designed_sequence(
    seq:     str,
    profile: list[dict],
    cols:    list[int] | None = None,
) -> list[dict]:
    """
    Score the residues of a designed sequence against the natural profile.

    `cols` restricts scoring to a list of 0-based alignment columns (see
    resolve_scoring_columns). Passing None scores every column, which
    reproduces the original behaviour.

    Metric A — entropy lookup:
        penalty_A = norm_H of that alignment column

        Careful: metric A is a property of the COLUMN and does not depend
        on which residue the design placed there. Across a design set that
        shares one gap pattern it is the same constant for every sequence
        and so carries no ranking information whatsoever. It is reported
        for inspection; rank on penalty_weighted instead.

    Metric B — distribution deviation:
        penalty_B = -log2( P(designed_AA | BLOSUM62-smoothed natural dist) )
                    divided by log2(20)

        Not bounded by 1. P is floored at 1e-8, so penalty_B can reach
        about 6.1 for a residue the family never tolerates.

    Weighted penalty — recommended for ranking:
        penalty_weighted = (1 - norm_H) * penalty_B

        The same substitution costs more at a conserved column than at a
        variable one. This is what metric A was reaching for; multiplying
        expresses it, whereas adding merely offsets every design equally.

    Combined = mean(A, B), retained only for comparability with earlier runs.
    """
    iter_cols = range(len(profile)) if cols is None else cols

    scores = []
    for col in iter_cols:
        if col >= len(seq) or col >= len(profile):
            continue
        prof = profile[col]
        aa   = seq[col]
        if aa in GAP_CHARS:
            continue

        aa     = aa.upper()
        norm_H = prof["norm_H"]

        p_aa          = max(prof["blosum_dist"].get(aa, 1e-8), 1e-8)
        neg_log_p     = -math.log2(p_aa)
        norm_neg_logp = neg_log_p / MAX_ENTROPY
        combined      = (norm_H + norm_neg_logp) / 2.0
        weighted      = (1.0 - norm_H) * norm_neg_logp

        scores.append({
            "col":                        col,
            "position":                   col + 1,
            "aa":                         aa,
            "gap_frac_natural":           prof["gap_frac"],
            "H_bits":                     prof["H_bits"],
            "penalty_A_entropy_lookup":   round(norm_H,          6),
            "penalty_B_distribution_dev": round(norm_neg_logp,   6),
            "penalty_weighted":           round(weighted,        6),
            "combined_penalty":           round(combined,        6),
            "p_aa_in_natural":            round(p_aa,            8),
        })
    return scores


def aggregate_scores(scores: list[dict]) -> dict:
    if not scores:
        return {}
    keys = [
        "penalty_A_entropy_lookup",
        "penalty_B_distribution_dev",
        "penalty_weighted",
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
    records:      list[tuple[str, str]],
    profile:      list[dict],
    prefix:       str,
    group:        str,
    cols:         list[int] | None = None,
    per_position: bool = False,
) -> list[dict]:
    """
    Score a group of sequences and write the result files.

    per_position controls the two heavyweight outputs (the full scores JSON
    and the long-format per-position CSV). They scale as sequences x columns
    and for a few thousand designs run to hundreds of megabytes, past what
    git will accept, so they are opt-in. The summary CSV is always written
    and is what you rank on.
    """
    group_slug = group.lower().strip().replace(" ", "_")
    all_results = []

    for header, seq in records:
        scores = score_designed_sequence(seq, profile, cols)
        agg    = aggregate_scores(scores)
        all_results.append({
            "header":       header,
            "sequence":     seq.replace("-", "").replace(".", ""),
            "aggregate":    agg,
            "per_position": scores,
        })

    # per-sequence JSON  (opt-in: large)
    if per_position:
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

    # long-format per-position CSV  (opt-in: large)
    if per_position:
        pos_csv = f"{prefix}_{group_slug}_per_position.csv"
        with open(pos_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([
                "sequence", "col", "position", "aa",
                "H_bits", "gap_frac_natural",
                "penalty_A_entropy_lookup",
                "penalty_B_distribution_dev",
                "penalty_weighted",
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
                        s["penalty_weighted"],
                        s["combined_penalty"],
                        s["p_aa_in_natural"],
                    ])
        print(f"  ✓ {group.title()} per-position CSV         → {pos_csv}")

    return all_results


def score_all_designed(
    records:      list[tuple[str, str]],
    profile:      list[dict],
    prefix:       str,
    cols:         list[int] | None = None,
    per_position: bool = False,
) -> list[dict]:
    return score_sequence_group(
        records, profile, prefix, "designed", cols, per_position
    )


def score_all_natural(
    records:      list[tuple[str, str]],
    profile:      list[dict],
    prefix:       str,
    cols:         list[int] | None = None,
    per_position: bool = False,
) -> list[dict]:
    return score_sequence_group(
        records, profile, prefix, "natural", cols, per_position
    )


# ══════════════════════════════════════════════════════════════════════════
# Calibration — turning a penalty into a percentile
# ══════════════════════════════════════════════════════════════════════════

def leave_one_out_calibration(
    records:      list[tuple[str, str]],
    cols:         list[int],
    max_gap_frac: float = 1.0,
    metric:       str = "mean_penalty_weighted",
) -> dict:
    """
    Score each natural sequence against a profile built from the others.

    Scoring the natural set against a profile built from all of it — which
    is what the original pipeline did — leaks: every sequence contributed
    to the distribution it is then measured against, so the scores come out
    flatteringly low and cannot serve as a reference. Holding each sequence
    out removes the leak and yields an honest reference distribution: the
    scores real members of this family earn.

    Column selection is kept fixed across the folds (from the full-set
    profile) so that every fold is scored on the same columns. Dropping one
    of 22 sequences does shift each column's gap fraction slightly, which
    would otherwise move the column set from fold to fold.
    """
    n = len(records)
    if n < 3:
        raise ValueError("Leave-one-out calibration needs at least 3 sequences.")

    folds = []
    for i, (header, seq) in enumerate(records):
        rest    = records[:i] + records[i + 1:]
        prof    = compute_natural_profile(rest, max_gap_frac=max_gap_frac)
        scores  = score_designed_sequence(seq, prof, cols)
        agg     = aggregate_scores(scores)
        folds.append({"header": header, **agg})

    vals = sorted(f[metric] for f in folds if metric in f)
    mean = sum(vals) / len(vals)
    var  = sum((v - mean) ** 2 for v in vals) / len(vals)
    return {
        "metric":       metric,
        "n_folds":      n,
        "scoring_cols": len(cols),
        "mean":         round(mean, 6),
        "sd":           round(var ** 0.5, 6),
        "min":          round(vals[0], 6),
        "max":          round(vals[-1], 6),
        "values":       [round(v, 6) for v in vals],
        "folds":        folds,
    }


def apply_calibration(
    all_results: list[dict],
    calib:       dict,
    metric:      str = "mean_penalty_weighted",
) -> None:
    """
    Annotate each scored sequence with its standing against the reference
    distribution, in place.

    percentile_vs_natural is the share of natural folds scoring BELOW this
    sequence, so 0 is better than every natural homolog and 100 is worse
    than all of them. With ~20 folds the resolution is about 5 points —
    read it as a band, not a number, and treat z_vs_natural as the finer
    ordering. worse_than_worst_natural is the blunt instrument: a design
    that trips it places a residue combination no real family member does.
    """
    ref = calib["values"]
    sd  = calib["sd"] or 1e-12
    mu  = calib["mean"]
    worst = calib["max"]

    for res in all_results:
        agg = res.get("aggregate") or {}
        val = agg.get(metric)
        if val is None:
            continue
        below = sum(1 for r in ref if r < val)
        agg["percentile_vs_natural"]    = round(100.0 * below / len(ref), 1)
        agg["z_vs_natural"]             = round((val - mu) / sd, 3)
        agg["worse_than_worst_natural"] = bool(val > worst)


def shuffle_null(
    records: list[tuple[str, str]],
    profile: list[dict],
    cols:    list[int],
    n_per:   int = 1,
    seed:    int = 0,
    metric:  str = "mean_penalty_weighted",
) -> dict:
    """
    Sanity control: permute each sequence's own residues among the scored
    columns and rescore.

    A shuffle preserves composition exactly and destroys only the pairing of
    residue to position, so it isolates whether the score responds to WHERE
    residues sit. If designs do not separate cleanly from their own shuffles,
    the metric is not measuring position-specific fit and nothing downstream
    of it means anything.
    """
    import random
    rng  = random.Random(seed)
    vals = []
    for _, seq in records:
        chars = [seq[c] for c in cols if c < len(seq)]
        for _ in range(n_per):
            shuffled = chars[:]
            rng.shuffle(shuffled)
            fake = dict(zip(cols, shuffled))
            scores = []
            for c in cols:
                if c not in fake:
                    continue
                aa = fake[c]
                if aa in GAP_CHARS:
                    continue
                prof   = profile[c]
                norm_H = prof["norm_H"]
                p_aa   = max(prof["blosum_dist"].get(aa.upper(), 1e-8), 1e-8)
                b      = -math.log2(p_aa) / MAX_ENTROPY
                scores.append((1.0 - norm_H) * b)
            if scores:
                vals.append(sum(scores) / len(scores))

    vals.sort()
    mean = sum(vals) / len(vals)
    var  = sum((v - mean) ** 2 for v in vals) / len(vals)
    return {
        "metric": metric,
        "n":      len(vals),
        "mean":   round(mean, 6),
        "sd":     round(var ** 0.5, 6),
        "min":    round(vals[0], 6),
        "max":    round(vals[-1], 6),
    }


# ══════════════════════════════════════════════════════════════════════════
# Template mapping — alignment columns ↔ structure residue numbers
# ══════════════════════════════════════════════════════════════════════════

class TemplateMap:
    """
    Two-way map between alignment columns and the template's residue numbers.

    Everything downstream of a resfile depends on this being right, and a
    wrong offset fails silently: Rosetta constrains the wrong residue at
    every position and reports no error. So the map is always built from a
    named sequence inside the alignment itself, never assumed.

    `start` is the residue number of the template's FIRST non-gap residue.
    Note that a PDB entry's author numbering need not match the order of
    its SEQRES records, and an RCSB FASTA gives the latter. Confirm `start`
    against the ATOM records of the PDB that Rosetta will read.
    """

    def __init__(self, header: str, aligned_seq: str, start: int = 1,
                 chain: str = "A"):
        self.header     = header
        self.aligned    = aligned_seq
        self.chain      = chain
        self.start      = start
        self.col_to_res: dict[int, int] = {}
        self.res_to_col: dict[int, int] = {}
        self.native:     dict[int, str] = {}

        resnum = start - 1
        for col, ch in enumerate(aligned_seq):
            if ch in GAP_CHARS:
                continue
            resnum += 1
            self.col_to_res[col] = resnum
            self.res_to_col[resnum] = col
            self.native[col] = ch.upper()

    @property
    def n_residues(self) -> int:
        return len(self.col_to_res)

    @property
    def res_range(self) -> tuple[int, int]:
        nums = sorted(self.res_to_col)
        return (nums[0], nums[-1]) if nums else (0, 0)

    def describe(self):
        lo, hi = self.res_range
        print(f"  Template            : {self.header[:60]}")
        print(f"  Template residues   : {self.n_residues}  "
              f"(numbered {lo}-{hi}, chain {self.chain})")


def find_template(
    records: list[tuple[str, str]],
    needle:  str,
) -> tuple[str, str]:
    """
    Locate the template sequence in the alignment by a header substring.

    Ambiguity is an error rather than a guess: picking the wrong one of two
    near-identical homologs would shift the whole numbering.
    """
    hits = [(h, s) for h, s in records if needle.lower() in h.lower()]
    if not hits:
        raise ValueError(
            f"--template {needle!r} matched no sequence header in the "
            f"alignment. Headers present:\n  " +
            "\n  ".join(h[:70] for h, _ in records[:8]) +
            ("\n  ..." if len(records) > 8 else "")
        )
    if len(hits) > 1:
        raise ValueError(
            f"--template {needle!r} matched {len(hits)} sequences; make it "
            f"specific:\n  " + "\n  ".join(h[:70] for h, _ in hits[:8])
        )
    return hits[0]


# ══════════════════════════════════════════════════════════════════════════
# Resfile I/O
# ══════════════════════════════════════════════════════════════════════════

RESFILE_COMMANDS = {
    "PIKAA", "NOTAA", "NATAA", "NATRO", "ALLAA", "ALLAAwc", "ALLAAxc",
    "POLAR", "APOLAR", "EX",
}


def looks_like_resfile(path: Path) -> bool:
    head = path.read_text(errors="replace")[:4000].upper()
    return "START" in head and any(c in head for c in ("PIKAA", "NATAA", "NATRO", "ALLAA"))


def parse_resfile(path: str | Path) -> dict[int, tuple[str, str, str]]:
    """
    Read a Rosetta resfile into {resnum: (chain, command, allowed_aas)}.

    Only the per-residue body after 'start' is read. The header block
    (default commands) is returned separately by read_resfile_default.
    """
    rows: dict[int, tuple[str, str, str]] = {}
    started = False
    for raw in Path(path).read_text(errors="replace").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if not started:
            if line.lower() == "start":
                started = True
            continue
        parts = line.split()
        if len(parts) < 3 or not parts[0].lstrip("-").isdigit():
            continue
        resnum, chain, cmd = int(parts[0]), parts[1], parts[2].upper()
        aas = parts[3].upper() if len(parts) > 3 else ""
        rows[resnum] = (chain, cmd, aas)
    return rows


def read_resfile_default(path: str | Path) -> str | None:
    """
    The default command that applies to residues NOT listed in the body.

    This is easy to get wrong and expensive when wrong: a resfile whose
    header reads '#NATRO' has that default COMMENTED OUT, so unlisted
    residues fall through to Rosetta's own default of designing them to
    any amino acid — the opposite of the intent.
    """
    for raw in Path(path).read_text(errors="replace").splitlines():
        stripped = raw.strip()
        if stripped.lower() == "start":
            return None
        line = stripped.split("#", 1)[0].strip()
        if line and line.split()[0].upper() in RESFILE_COMMANDS:
            return line
    return None


# ══════════════════════════════════════════════════════════════════════════
# Resfile generation from the entropy profile
# ══════════════════════════════════════════════════════════════════════════

def nucleus_alphabet(
    dist:   dict[str, float],
    top_p:  float,
    min_aa: int,
    max_aa: int,
    native: str | None = None,
) -> tuple[list[str], float]:
    """
    The smallest set of amino acids whose natural probability sums to top_p.

    A fixed probability threshold gives erratic set sizes — one residue at a
    conserved column, fifteen at a variable one. Taking a cumulative mass
    instead lets the set size follow the column's own tolerance, which is
    the quantity of interest.

    The template's native residue is always included: excluding it would let
    Rosetta be forced off a wild-type residue that is known to work.
    """
    ranked = sorted(AA_ORDER, key=lambda a: -dist.get(a, 0.0))
    chosen, cum = [], 0.0
    for aa in ranked:
        chosen.append(aa)
        cum += dist.get(aa, 0.0)
        if cum >= top_p and len(chosen) >= min_aa:
            break
    if len(chosen) > max_aa:
        chosen = chosen[:max_aa]
        cum = sum(dist.get(a, 0.0) for a in chosen)
    if native and native in AA_SET and native not in chosen:
        chosen.append(native)
    chosen.sort(key=lambda a: -dist.get(a, 0.0))
    return chosen, cum


def build_resfile(
    profile:        list[dict],
    tmap:           TemplateMap,
    conserved_max_h: float = 0.15,
    top_p:          float = 0.90,
    min_aa:         int = 1,
    max_aa:         int = 20,
    lock:           set[int] | None = None,
    lock_cys:       bool = True,
) -> list[dict]:
    """
    Turn the entropy profile into a per-residue design specification.

    Entropy chooses the MODE — whether a position is designed at all.
    The distribution chooses the ALPHABET — what is permitted where it is.
    """
    lock = lock or set()
    decisions = []

    for col, prof in enumerate(profile):
        if col not in tmap.col_to_res:
            continue                       # template has a gap here
        resnum = tmap.col_to_res[col]
        native = tmap.native[col]
        norm_H = prof["norm_H"]

        if resnum in lock:
            mode, aas, reason = "NATAA", "", "locked by --lock"
        elif lock_cys and native == "C":
            mode, aas, reason = "NATAA", "", "cysteine (possible disulfide)"
        elif not prof.get("scorable", True):
            mode, aas, reason = "NATAA", "", f"gappy column ({prof['gap_frac']:.2f})"
        elif norm_H <= conserved_max_h:
            mode, aas, reason = "NATAA", "", f"conserved (norm_H {norm_H:.3f})"
        else:
            alpha, cum = nucleus_alphabet(
                prof["blosum_dist"], top_p, min_aa, max_aa, native
            )
            if len(alpha) <= 1:
                mode, aas, reason = "NATAA", "", "single residue tolerated"
            else:
                mode, aas = "PIKAA", "".join(alpha)
                reason = f"top_p {cum:.2f}, {len(alpha)} aa"

        decisions.append({
            "resnum":   resnum,
            "col":      col,
            "position": col + 1,
            "chain":    tmap.chain,
            "native":   native,
            "norm_H":   round(norm_H, 6),
            "gap_frac": prof["gap_frac"],
            "command":  mode,
            "allowed":  aas,
            "n_allowed": len(aas) if aas else 1,
            "reason":   reason,
        })
    return decisions


def write_resfile(decisions: list[dict], path: str, default_cmd: str = "NATRO") -> None:
    with open(path, "w") as f:
        f.write("# generated by entropy_pipeline.py make-resfile\n")
        f.write("# entropy chooses the mode, the natural distribution the alphabet\n")
        f.write("# verify residue numbering against the PDB that Rosetta will read\n")
        f.write(f"{default_cmd}\n")
        f.write("start\n")
        for d in decisions:
            if d["command"] == "PIKAA":
                f.write(f"{d['resnum']:<4d} {d['chain']}  PIKAA  {d['allowed']}\n")
            else:
                f.write(f"{d['resnum']:<4d} {d['chain']}  {d['command']}\n")
    print(f"  ✓ Resfile → {path}")


def write_resfile_report(decisions: list[dict], path: str) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(decisions[0].keys()))
        w.writeheader()
        w.writerows(decisions)
    print(f"  ✓ Per-residue rationale → {path}")


def compare_resfiles(decisions: list[dict], reference: str) -> None:
    """
    Diff a generated specification against a hand-built one.

    Agreement is evidence the generator works; disagreement says either the
    generator or the manual choice is wrong, and both are worth knowing.
    Containment is the more informative direction: if the manual alphabet
    sits inside the generated one, the generator is permissive enough to
    have allowed every choice that was actually made by hand.
    """
    ref = parse_resfile(reference)
    gen = {d["resnum"]: (d["command"], d["allowed"]) for d in decisions}
    dflt = read_resfile_default(reference)

    print("\n  ───────────────────────────────────────────────────────")
    print(f"  Comparison against {reference}")
    print("  ───────────────────────────────────────────────────────")
    print(f"  Reference positions : {len(ref)}")
    print(f"  Generated positions : {len(gen)}")

    if dflt is None:
        print("  !! The reference sets no default command before 'start', so")
        print("     residues it does not list fall through to Rosetta's own")
        print("     default and are designable to any amino acid. If the header")
        print("     has a commented-out '#NATRO', that is likely unintended.")
    else:
        print(f"  Reference default   : {dflt}")

    both  = sorted(set(ref) & set(gen))
    print(f"  In both             : {len(both)}")
    print(f"  Reference only      : {len(set(ref) - set(gen))}")
    print(f"  Generated only      : {len(set(gen) - set(ref))}")
    if not both:
        print("  No shared positions — the numbering almost certainly differs.")
        print("  ───────────────────────────────────────────────────────")
        return

    same_mode = sum(1 for r in both if ref[r][1] == gen[r][0])
    print(f"  Same command        : {same_mode}/{len(both)} "
          f"({100*same_mode/len(both):.0f}%)")

    pairs = [(set(ref[r][2]), set(gen[r][1])) for r in both
             if ref[r][1] == "PIKAA" and gen[r][0] == "PIKAA"]
    if pairs:
        jac = [len(a & b) / len(a | b) for a, b in pairs if a | b]
        print(f"  PIKAA on both sides : {len(pairs)}")
        print(f"    mean overlap      : {sum(jac)/len(jac):.2f} (Jaccard)")
        print(f"    manual ⊆ generated: {sum(1 for a,b in pairs if a <= b)}"
              f"/{len(pairs)}")
        print(f"    generated ⊆ manual: {sum(1 for a,b in pairs if b <= a)}"
              f"/{len(pairs)}")
        ref_sz = sum(len(a) for a, _ in pairs) / len(pairs)
        gen_sz = sum(len(b) for _, b in pairs) / len(pairs)
        print(f"    mean alphabet size: manual {ref_sz:.1f}  generated {gen_sz:.1f}")
    print("  ───────────────────────────────────────────────────────")


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

def _report_calibration(calib: dict, shuf: dict | None = None):
    print("\n  ───────────────────────────────────────────────────────")
    print("  Leave-one-out reference distribution")
    print("  ───────────────────────────────────────────────────────")
    print(f"  Metric              : {calib['metric']}")
    print(f"  Folds (held out)    : {calib['n_folds']}")
    print(f"  Scoring columns     : {calib['scoring_cols']}")
    print(f"  Natural mean ± sd   : {calib['mean']:.4f} ± {calib['sd']:.4f}")
    print(f"  Natural range       : {calib['min']:.4f} – {calib['max']:.4f}")
    if shuf:
        sep = (shuf["mean"] - calib["mean"]) / (calib["sd"] or 1e-12)
        print(f"  Shuffle-null mean   : {shuf['mean']:.4f} ± {shuf['sd']:.4f}")
        print(f"  Separation          : {sep:+.1f} natural sd")
        if sep < 1.0:
            print("  !! Shuffles score about as well as the real sequences.")
            print("     The metric is not responding to residue placement —")
            print("     fix that before trusting any ranking below.")
    print("  ───────────────────────────────────────────────────────")


def _build_tmap(records, args):
    """Build a TemplateMap if --template was given, else None."""
    needle = getattr(args, "template", None)
    if not needle:
        return None
    header, seq = find_template(records, needle)
    tmap = TemplateMap(
        header, seq,
        start=getattr(args, "template_start", 1),
        chain=getattr(args, "chain", "A"),
    )
    tmap.describe()
    return tmap


def _resolve_cols_for_designs(profile, designed, args, tmap=None):
    """Work out which columns to score, honouring --positions."""
    spec = getattr(args, "positions", None)
    if spec is not None and spec.strip().lower() == "auto":
        varied = detect_variable_columns(designed)
        print(f"\n  --positions auto: {len(varied)} of {len(profile)} columns "
              f"differ between designs")
        if not varied:
            print("    none vary — falling back to every column")
            positions = None
        else:
            positions = varied
    else:
        positions = parse_positions(spec, len(profile), tmap)
    return resolve_scoring_columns(profile, positions)


def cmd_natural(args):
    print(f"\n[Stage 1] Natural sequences: {args.fasta}")
    records = read_fasta(args.fasta)
    print(f"  Loaded {len(records)} sequences  |  alignment length {len(records[0][1])}")
    profile = compute_natural_profile(records, max_gap_frac=args.max_gap_frac)
    n_drop = sum(1 for c in profile if not c["scorable"])
    if n_drop:
        print(f"  Columns above --max-gap-frac {args.max_gap_frac}: "
              f"{n_drop} marked unscorable")
    print_natural_summary(profile, len(records))
    save_natural_profile(profile, args.out)

    if args.loo:
        positions = parse_positions(getattr(args, "positions", None), len(profile))
        cols  = resolve_scoring_columns(profile, positions)
        calib = leave_one_out_calibration(
            records, cols, max_gap_frac=args.max_gap_frac
        )
        shuf  = shuffle_null(records, profile, cols, n_per=20)
        _report_calibration(calib, shuf)
        path = f"{args.out}_calibration.json"
        with open(path, "w") as f:
            json.dump({"leave_one_out": calib, "shuffle_null": shuf}, f, indent=2)
        print(f"  ✓ Calibration → {path}")

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

    tmap = _build_tmap(all_records, args)
    cols = _resolve_cols_for_designs(profile, designed, args, tmap)

    all_results = score_all_designed(
        designed, profile, args.out, cols, args.per_position
    )
    natural_results = score_all_natural(
        natural, profile, args.out, cols, args.per_position
    )
    print_designed_summary(all_results)
    print_natural_penalty_summary(natural_results)

    if args.calibration:
        with open(args.calibration) as f:
            blob = json.load(f)
        calib = blob.get("leave_one_out", blob)

        # A percentile only means something if the reference distribution was
        # built on the same columns the designs were scored on. In the
        # two-step workflow it is easy to calibrate over the whole alignment
        # and then score designs on a subset, which silently compares two
        # different quantities.
        calib_cols = calib.get("scoring_cols")
        if calib_cols is not None and calib_cols != len(cols):
            print(f"\n  !! Calibration was built on {calib_cols} columns but these")
            print(f"     designs were scored on {len(cols)}. Percentiles and")
            print( "     z-scores would compare different quantities, so they are")
            print( "     being skipped. Re-run 'natural --loo' with the same")
            print( "     --positions, or use the 'both --loo' subcommand, which")
            print( "     keeps the two consistent by construction.")
            args.calibration = None

    if args.calibration:
        apply_calibration(all_results, calib)
        apply_calibration(natural_results, calib)
        _report_calibration(calib, blob.get("shuffle_null"))
        flagged = sum(
            1 for r in all_results
            if r["aggregate"].get("worse_than_worst_natural")
        )
        print(f"\n  Designs worse than every natural homolog: "
              f"{flagged} / {len(all_results)}")
        # rewrite the summary CSV so the calibrated fields land in it
        rows = [{"sequence": r["header"], **r["aggregate"]} for r in all_results]
        if rows:
            sum_csv = f"{args.out}_designed_summary.csv"
            with open(sum_csv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
            print(f"  ✓ Calibrated summary → {sum_csv}")

    if args.plot:
        plot_designed_scores(all_results, profile, args.out)


def cmd_both(args):
    # Stage 1
    print(f"\n[Stage 1] Natural sequences: {args.natural}")
    nat_records = read_fasta(args.natural)
    print(f"  Loaded {len(nat_records)} sequences  |  alignment length {len(nat_records[0][1])}")
    profile = compute_natural_profile(nat_records, max_gap_frac=args.max_gap_frac)
    n_drop = sum(1 for c in profile if not c["scorable"])
    if n_drop:
        print(f"  Columns above --max-gap-frac {args.max_gap_frac}: "
              f"{n_drop} marked unscorable")
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

    tmap = _build_tmap(nat_records, args)
    cols = _resolve_cols_for_designs(profile, designed, args, tmap)

    all_results = score_all_designed(
        designed, profile, args.out, cols, args.per_position
    )
    natural_results = score_all_natural(
        natural, profile, args.out, cols, args.per_position
    )
    print_designed_summary(all_results)
    print_natural_penalty_summary(natural_results)

    if args.loo:
        calib = leave_one_out_calibration(
            nat_records, cols, max_gap_frac=args.max_gap_frac
        )
        shuf = shuffle_null(nat_records, profile, cols, n_per=20)
        apply_calibration(all_results, calib)
        apply_calibration(natural_results, calib)
        _report_calibration(calib, shuf)
        flagged = sum(
            1 for r in all_results
            if r["aggregate"].get("worse_than_worst_natural")
        )
        print(f"\n  Designs worse than every natural homolog: "
              f"{flagged} / {len(all_results)}")
        path = f"{args.out}_calibration.json"
        with open(path, "w") as f:
            json.dump({"leave_one_out": calib, "shuffle_null": shuf}, f, indent=2)
        print(f"  ✓ Calibration → {path}")
        rows = [{"sequence": r["header"], **r["aggregate"]} for r in all_results]
        if rows:
            sum_csv = f"{args.out}_designed_summary.csv"
            with open(sum_csv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
            print(f"  ✓ Calibrated summary → {sum_csv}")

    if args.plot:
        plot_designed_scores(all_results, profile, args.out)


def cmd_make_resfile(args):
    print(f"\n[make-resfile] Natural alignment: {args.fasta}")
    records = read_fasta(args.fasta)
    print(f"  Loaded {len(records)} sequences  |  alignment length {len(records[0][1])}")

    header, seq = find_template(records, args.template)
    tmap = TemplateMap(header, seq, start=args.template_start, chain=args.chain)
    tmap.describe()

    if args.profile:
        print(f"  Reusing profile     : {args.profile}")
        with open(args.profile) as f:
            profile = json.load(f)
        if len(profile) != len(records[0][1]):
            raise ValueError(
                f"Profile has {len(profile)} columns but the alignment has "
                f"{len(records[0][1])}. They must come from the same alignment."
            )
    else:
        profile = compute_natural_profile(records, max_gap_frac=args.max_gap_frac)

    lock = set()
    if args.lock:
        lock = {int(x) for x in args.lock.replace(" ", "").split(",") if x}
        outside = sorted(r for r in lock if r not in tmap.res_to_col)
        if outside:
            lo, hi = tmap.res_range
            print(f"  !! --lock lists residues outside the template "
                  f"({lo}-{hi}): {outside}")

    decisions = build_resfile(
        profile, tmap,
        conserved_max_h=args.conserved_max_h,
        top_p=args.top_p,
        min_aa=args.min_aa,
        max_aa=args.max_aa,
        lock=lock,
        lock_cys=not args.no_lock_cys,
    )

    n_design = sum(1 for d in decisions if d["command"] == "PIKAA")
    n_fixed  = len(decisions) - n_design
    sizes    = [d["n_allowed"] for d in decisions if d["command"] == "PIKAA"]

    print("\n  ───────────────────────────────────────────────────────")
    print("  Design space")
    print("  ───────────────────────────────────────────────────────")
    print(f"  Template residues   : {len(decisions)}")
    print(f"  Designable (PIKAA)  : {n_design}")
    print(f"  Held fixed (NATAA)  : {n_fixed}")
    if sizes:
        print(f"  Alphabet size       : mean {sum(sizes)/len(sizes):.1f}  "
              f"min {min(sizes)}  max {max(sizes)}")
        log10 = sum(math.log10(s) for s in sizes)
        print(f"  Sequence space      : ~1e{log10:.0f}")
    reasons = Counter(d["reason"].split("(")[0].split(",")[0].strip()
                      for d in decisions if d["command"] != "PIKAA")
    for why, n in reasons.most_common():
        print(f"    fixed: {why:<32s} {n}")
    print("  ───────────────────────────────────────────────────────")

    out = args.out if args.out.endswith(".resfile") else f"{args.out}.resfile"
    write_resfile(decisions, out, default_cmd=args.default_command)
    write_resfile_report(decisions, out.replace(".resfile", "_rationale.csv"))

    if args.compare:
        compare_resfiles(decisions, args.compare)

    print("\n  Before running Rosetta: confirm the residue numbering above")
    print("  matches the ATOM records of the PDB you will pass it.")


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
    p1.add_argument("--max-gap-frac", type=float, default=0.3,
        help="Mark columns with a natural gap fraction above this as\n"
             "unscorable. A column where most homologs are absent has its\n"
             "distribution estimated from a few residues and produces\n"
             "unreliable extremes.  (default: 0.3)")
    p1.add_argument("--positions", default=None,
        help="Alignment columns to score:\n"
             "  auto          only columns that differ between designs\n"
             "                (the designable positions)\n"
             "  all           every column  (default)\n"
             "  80,94-99      explicit 1-based alignment positions\n"
             "  cols.txt      file, one position per line\n"
             "These are alignment columns, not PDB residue numbers.")
    p1.add_argument("--loo", action="store_true",
        help="Leave-one-out calibration: score each natural sequence against\n"
             "a profile built from the others, plus a shuffle control.\n"
             "Writes <out>_calibration.json for use with 'designed'.")

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
    p2.add_argument("--positions", default=None,
        help="Alignment columns to score:\n"
             "  auto          only columns that differ between designs\n"
             "                (the designable positions)\n"
             "  all           every column  (default)\n"
             "  80,94-99      explicit 1-based alignment positions\n"
             "  cols.txt      file, one position per line\n"
             "These are alignment columns, not PDB residue numbers.")
    p2.add_argument("--calibration", default=None,
        help="Calibration JSON from 'natural --loo'. Adds percentile,\n"
             "z-score and worse-than-worst-natural flags to the summary.")
    p2.add_argument("--per-position", action="store_true",
        help="Also write the full per-residue scores JSON and long-format\n"
             "CSV. These scale as sequences x columns and reach hundreds of\n"
             "megabytes for a few thousand designs, so they are off by default.")
    p2.add_argument("--template", default=None,
        help="Substring of the template sequence's header in the alignment\n"
             "(e.g. 'P00590'). Required to use a resfile with --positions.")
    p2.add_argument("--template-start", type=int, default=1,
        help="Residue number of the template's first non-gap residue\n"
             "(default: 1). Must match the PDB numbering Rosetta will read.")


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
    p3.add_argument("--max-gap-frac", type=float, default=0.3,
        help="Mark columns with a natural gap fraction above this as\n"
             "unscorable. A column where most homologs are absent has its\n"
             "distribution estimated from a few residues and produces\n"
             "unreliable extremes.  (default: 0.3)")
    p3.add_argument("--positions", default=None,
        help="Alignment columns to score:\n"
             "  auto          only columns that differ between designs\n"
             "                (the designable positions)\n"
             "  all           every column  (default)\n"
             "  80,94-99      explicit 1-based alignment positions\n"
             "  cols.txt      file, one position per line\n"
             "These are alignment columns, not PDB residue numbers.")
    p3.add_argument("--loo", action="store_true",
        help="Leave-one-out calibration of the natural set, and annotate\n"
             "every design with its standing against that distribution.")
    p3.add_argument("--per-position", action="store_true",
        help="Also write the large per-residue JSON and CSV outputs.")
    p3.add_argument("--template", default=None,
        help="Substring of the template sequence's header in the alignment\n"
             "(e.g. 'P00590'). Required to use a resfile with --positions.")
    p3.add_argument("--template-start", type=int, default=1,
        help="Residue number of the template's first non-gap residue\n"
             "(default: 1). Must match the PDB numbering Rosetta will read.")


    # ── make-resfile ──
    p4 = sub.add_parser("make-resfile",
        help="Turn the entropy profile into a Rosetta resfile")
    p4.add_argument("--fasta", required=True,
        help="MAFFT-aligned natural sequences, INCLUDING the template")
    p4.add_argument("--template", required=True,
        help="Substring of the template sequence's header in the alignment\n"
             "(e.g. 'P00590'). The template must be in the alignment so the\n"
             "column-to-residue map is derived rather than assumed.")
    p4.add_argument("--template-start", type=int, default=1,
        help="Residue number of the template's first non-gap residue\n"
             "(default: 1). Confirm against the ATOM records of the PDB\n"
             "Rosetta will read — an RCSB FASTA gives SEQRES order, which\n"
             "need not match author numbering.")
    p4.add_argument("--chain", default="A",
        help="Chain identifier written into the resfile  (default: A)")
    p4.add_argument("--profile", default=None,
        help="Reuse an existing *_profile.json instead of recomputing.\n"
             "Must come from the same alignment.")
    p4.add_argument("--out", default="designed",
        help="Output path; '.resfile' is appended if absent  (default: designed)")
    p4.add_argument("--max-gap-frac", type=float, default=0.3,
        help="Columns gappier than this are held fixed  (default: 0.3)")
    p4.add_argument("--conserved-max-h", type=float, default=0.15,
        help="Hold a position fixed when its normalised entropy is at or\n"
             "below this. Raise it to design more conservatively, lower it\n"
             "to open up more positions.  (default: 0.15)")
    p4.add_argument("--top-p", type=float, default=0.90,
        help="Allow the smallest set of amino acids whose natural\n"
             "probability sums to this  (default: 0.90)")
    p4.add_argument("--min-aa", type=int, default=1,
        help="Minimum alphabet size at a designed position  (default: 1)")
    p4.add_argument("--max-aa", type=int, default=20,
        help="Maximum alphabet size at a designed position  (default: 20)")
    p4.add_argument("--lock", default=None,
        help="Comma-separated residue numbers to hold fixed regardless of\n"
             "entropy — the catalytic triad, anything mechanistically\n"
             "essential. In the template's numbering, so verify it against\n"
             "the structure rather than taking numbers from a paper.")
    p4.add_argument("--no-lock-cys", action="store_true",
        help="Allow design at cysteines. Off by default, since mutating one\n"
             "half of a disulfide is rarely intended.")
    p4.add_argument("--default-command", default="NATRO",
        help="Default command for residues not listed  (default: NATRO).\n"
             "Written uncommented: a commented-out default silently lets\n"
             "Rosetta design every unlisted residue.")
    p4.add_argument("--compare", default=None,
        help="An existing resfile to diff the generated one against")

    args = parser.parse_args()
    {
        "natural":      cmd_natural,
        "designed":     cmd_designed,
        "both":         cmd_both,
        "make-resfile": cmd_make_resfile,
    }[args.cmd](args)
    print("Done.\n")


if __name__ == "__main__":
    main()
