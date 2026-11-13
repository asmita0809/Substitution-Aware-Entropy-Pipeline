# Substitution-Aware Entropy Pipeline

Scores designed protein sequences against the evolutionary tolerance of their
natural homologs, using BLOSUM62-weighted per-position Shannon entropy.

Built for a cutinase design campaign — cutinase (PDB **1CEX**, *Fusarium solani*)
is an ester hydrolase of interest for polyester/PET degradation — but the
pipeline is sequence-family agnostic. Point it at any MAFFT alignment.

## The idea

A protein design pipeline will happily hand you thousands of candidate
sequences with no indication of which ones violate constraints that evolution
has been enforcing for millions of years. This tool asks a narrower, cheaper
question than a folding or stability calculation:

> At each position, how much variation does the natural family actually
> tolerate — and does this designed sequence respect that?

Positions where natural homologs are invariant are positions where a
substitution is likely to be costly. Positions that drift freely across the
family are positions where a design can be bold for free.

### Why BLOSUM62-weighted, not raw frequency

Plain per-column frequency counting treats every substitution as equally
foreign. With 22 natural sequences, a column containing only `I` will score a
designed `V` exactly as harshly as a designed `D` — even though I→V is a
conservative aliphatic swap and I→D is not.

So each observed residue is first spread over all 20 amino acids via a softmax
over its BLOSUM62 row, and the column distribution is the mean of those spread
distributions. Shannon entropy is then taken over that smoothed distribution.
The effect is twofold: chemically similar substitutions are pre-credited, and
sparse columns stop producing false zero-probability cliffs.

## The two metrics

Each residue of each designed sequence gets scored on two independent axes,
then averaged:

| Metric | Definition | What it catches |
|---|---|---|
| **A — column entropy lookup** | `norm_H` at that alignment column, in [0,1] | Mutating a position the family holds rigid |
| **B — distribution deviation** | `-log2 P(designed AA \| BLOSUM natural dist)` | Choosing a residue the family specifically avoids |
| **Combined** | `mean(A, B)` | — |

The two can disagree, and that disagreement is informative. A high-entropy
column with a chemically bizarre residue scores low on A and high on B. A
rigidly conserved column where the design picked the conserved residue anyway
scores low on both. Metric A is a property of the *position*; metric B is a
property of the *choice made at it*.

Lower is better throughout. Scores are penalties, not fitness.

## Install

```bash
git clone <your-repo-url>
cd entropy-pipeline
```

The core pipeline runs on the **Python standard library alone** (3.9+) — no
install step. `matplotlib` and `numpy` are imported lazily, only when you pass
`--plot`:

```bash
pip install -r requirements.txt   # only needed for figures / the notebook
```

[MAFFT](https://mafft.cbrc.jp/alignment/software/) is required to produce the
alignments, but is not called by this script — you align first, then feed the
result in.

## Quickstart

A 22-sequence natural cutinase alignment ships in `data/`, so Stage 1 runs
with no setup:

```bash
python entropy_pipeline.py natural \
    --fasta data/natural_cutinase_msa_22.fasta \
    --out natural --plot
```

That writes `natural_profile.json` (the reusable profile), `natural_profile.csv`
(per-position frequencies, BLOSUM distributions and entropy), and a PNG.

## Full workflow

```bash
# 1 — align the natural homologs
mafft --auto natural.fasta > natural_aligned.fasta

# 2 — add designed sequences into that existing alignment
#     (--add keeps the natural columns fixed, which is what makes
#      column indices comparable between the two stages)
mafft --add designed.fasta --reorder natural_aligned.fasta > combined.fasta

# 3 — build the natural entropy profile
python entropy_pipeline.py natural --fasta natural_aligned.fasta --out natural --plot

# 4 — score the designs
python entropy_pipeline.py designed \
    --fasta combined.fasta \
    --profile natural_profile.json \
    --design-prefix "design_" \
    --out results --plot
```

Or collapse steps 3–4:

```bash
python entropy_pipeline.py both \
    --natural natural_aligned.fasta \
    --combined combined.fasta \
    --design-prefix "design_" \
    --out results --plot
```

`--design-prefix` is matched case-insensitively against the full FASTA header,
which is how natural and designed records are separated out of the single
combined file. Sanity-check the split report it prints before trusting the
scores — a prefix that matches nothing silently yields zero designs.

## Outputs

| File | Contents |
|---|---|
| `*_profile.json` | Per-column BLOSUM distribution + entropy. Reusable across runs. |
| `*_profile.csv` | Same, flat: `position, gap_frac, n_seqs, H_bits, norm_H, freq_*, blosum_*` |
| `*_scores.json` | Full per-residue scores for every sequence |
| `*_per_position.csv` | Long-format per-residue table |
| `*_summary.csv` | One row per sequence: mean/sum/max/min of each metric |

**Scale warning.** `*_scores.json` and `*_per_position.csv` grow with
sequences × alignment columns. The 5,150-design run behind this repo produced a
319 MB JSON and an 87 MB CSV — both over or near GitHub's 100 MB file ceiling,
and both excluded by `.gitignore`. Work from `*_summary.csv` for ranking and
reach for the per-position files only when inspecting specific candidates.

## Notebook

`entropy_pipeline.ipynb` is a cell-by-cell version of the same code, plus a
`logomaker` sequence-logo stage that the CLI does not have. Outputs are
cleared — run top to bottom. The config cell near the top sets the input
filenames and design prefix in one place.

## Repository layout

```
entropy_pipeline.py        # CLI: natural / designed / both
entropy_pipeline.ipynb     # notebook twin + sequence logos
data/
  natural_cutinase_msa_22.fasta   # aligned, 22 homologs — Stage 1 demo input
  natural_cutinase_23.fasta       # 23-sequence set
  fusarium_msa_27.fasta           # 27-sequence Fusarium alignment
  1cex_reference.fasta            # PDB 1CEX reference sequence
examples/
  profiles/                # precomputed Stage 1 output on the 22-seq MSA
  figures/                 # entropy profile + sequence logos
```

## Example result

On the bundled 22-sequence cutinase alignment: **276 alignment columns**, mean
normalized entropy **0.393**, with **60 columns below 0.2** — roughly a fifth of
the alignment is strongly conserved and should be treated as near-untouchable
by a design protocol. See `examples/figures/natural_entropy_profile.png`.

## Limitations

- Entropy is computed per column, independently. Coupled positions — salt
  bridges, packing pairs, catalytic triads acting together — are invisible to
  it. A design can score well here and still be broken.
- Results depend on the natural set. 22 sequences is a thin sample; conservation
  at a column may reflect shared ancestry rather than functional constraint.
  Broader, less redundant homolog sets give better-calibrated entropy.
- No structural term. Burial, secondary structure and solvent accessibility all
  modulate how costly a substitution actually is, and none are used.
- Gaps are skipped rather than modeled, so indel-heavy designs are scored on
  fewer positions than their length suggests.

## Background reading

Not included in this repository (copyrighted); listed for provenance of the
approach:

- Chan et al. (2004) — Relationship between local structural entropy and protein thermostability
- Bae et al. (2008) — Bioinformatic method for protein thermal stabilization by structural entropy optimization
- Borgo & Havranek (2012) — Automated selection of stabilizing mutations in designed and natural proteins
- Isakova et al. (2026) — Descent from a common ancestor restricts exploration of protein sequence space
- Additional reading on relative-entropy protein design, evolution-entropy hybrid residue ranking, local entropy in proteins, and thermal stability in plastic-degrading enzymes

## Status

Exploratory research code from a working design campaign, cleaned up and made
reproducible. It is not a validated predictor — treat the scores as a triage
filter for prioritizing candidates, not as a stability or activity prediction.

## License

MIT — see `LICENSE`.
