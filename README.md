# Substitution-Aware Entropy Pipeline

Scores designed protein sequences against the evolutionary tolerance of their
natural homologs, using BLOSUM62-weighted per-position Shannon entropy.

Built for a cutinase design campaign — cutinase (PDB **1CEX**, *Fusarium solani*)
is an ester hydrolase of interest for polyester/PET degradation — but the
pipeline is sequence-family agnostic. Point it at any MAFFT alignment.

## The idea

A design pipeline will happily hand you thousands of candidate sequences with
no indication of which ones violate constraints evolution has been enforcing
for millions of years. This tool asks a narrower, cheaper question than a
folding or stability calculation:

> At each position, how much variation does the natural family actually
> tolerate — and does this designed sequence respect that?

Positions where natural homologs are invariant are positions where a
substitution is likely to be costly. Positions that drift freely across the
family are positions where a design can be bold for free.

### Why BLOSUM62-weighted, not raw frequency

Plain per-column frequency counting treats every substitution as equally
foreign. With 22 natural sequences, a column containing only `I` scores a
designed `V` exactly as harshly as a designed `D` — though I→V is a
conservative aliphatic swap and I→D is not.

So each observed residue is first spread over all 20 amino acids via a softmax
over its BLOSUM62 row, and the column distribution is the mean of those spread
distributions. Shannon entropy is taken over that smoothed distribution.
Chemically similar substitutions are pre-credited, and sparse columns stop
producing false zero-probability cliffs.

## Relation to existing methods

**This is a lightweight reimplementation of an established idea, not a new
method.** Restricting protein design using conservation derived from an
alignment of homologs is well developed, implemented in widely used servers,
and experimentally validated at a scale nothing here approaches. If you need
design output you can act on, use those tools. If you want to see and tune
every step of the reasoning, use this one.

**PROSS** (Goldenzweig et al., *Mol Cell* 2016; PROSS 2 server, Weinstein et
al., *Bioinformatics* 2021) builds a PSSM from homologs and combines it with
Rosetta all-atom energies; the current server accepts an uploaded MSA directly.
Its community evaluation is unusually thorough — twelve independent
laboratories ran it on 14 unrelated targets without developer support, and nine
of ten tested targets gained thermal stability (Peleg et al., *JMB* 2021).

**FuncLib** (Khersonsky et al., *Mol Cell* 2018) is the closest relative to
`make-resfile`. It constrains sequence space position by position using
conservation (PSSM) together with ΔΔG mutational scanning, producing a
per-position list of tolerated residues, then enumerates combinations with
Rosetta. That per-position tolerated alphabet is the same object `make-resfile`
emits — the difference is that FuncLib adds a structural energy term and this
does not.

**htFuncLib** (Weinstein et al., *Nat Commun* 2023) extends that with machine
learning to design spaces whose mutations combine without incompatible
interactions, recovering thousands of functional GFP variants. That directly
addresses the per-column independence listed under Known Limitations below.

Using entropy specifically, rather than a PSSM cutoff, also has its own
lineage — Chan et al. (2004), Bae et al. (2008) and Borgo & Havranek (2012),
cited at the end.

And for *ranking* designs, learned models have overtaken profile methods
outright. Inverse-folding models such as ProteinMPNN and ESM-IF substantially
outperform sequence-statistics baselines in experimental benchmarks. A
BLOSUM-smoothed profile over 22 sequences will not compete with them.

### What this adds

Nothing algorithmic. What it does offer:

- **Calibration surfaced to the user.** The leave-one-out reference
  distribution and the shuffle control tell you what scores real family
  members earn, and whether the metric responds to residue placement at all.
  Tools in this space generally don't hand you that.
- **Every decision carries its reason.** `*_rationale.csv` records why each
  residue was opened or held fixed, which is how the residue-160 case below
  was found.
- **No dependencies, one second.** The core runs on the Python standard
  library, so parameters can be swept locally rather than resubmitted to a
  server.
- **Numbering checked rather than assumed**, with the offset verified against
  a named template sequence inside the alignment.

If you use this in published work, cite the methods above — not this
repository.

## Scoring

Each scored residue gets three numbers:

| Metric | Definition | Range | Use |
|---|---|---|---|
| **A — column entropy** | `norm_H` at that column | [0, 1] | diagnostic only |
| **B — distribution deviation** | `−log2 P(aa \| BLOSUM dist) / log2(20)` | [0, ~6.1] | component |
| **Weighted** | `(1 − norm_H) · B` | — | **rank on this** |
| *Combined (legacy)* | `mean(A, B)` | — | comparability with old runs |

**Metric A cannot rank designs, and the legacy combined score inherits that.**
A is a property of the *column*, not of the residue placed in it. Across a
design set that shares one gap pattern, every sequence averages A over the same
columns and receives the same constant. On the 5,150-design run behind this
repo, metric A was identical to six decimal places for every design, and
`mean(A, B)` therefore spanned only 4.2% of its own mean — noise.

The weighted penalty expresses what A was reaching for: the same substitution
costs more at a conserved column than at a variable one. Multiplying encodes
that; adding merely offsets every design equally.

Lower is better throughout. These are penalties, not fitness.

## Column selection

Two filters decide which columns are scored, and both matter more than the
choice of metric.

**`--max-gap-frac` (default 0.3)** marks columns where too many natural
sequences carry a gap as unscorable. Such a column's distribution rests on a
handful of residues and will cheerfully emit the largest penalties in a run out
of pure sampling noise. In the original cutinase run the three worst-scoring
variable positions were columns that are 73%, 55% and 82% gaps in the natural
set — estimated from 6, 10 and 4 sequences respectively.

**`--positions auto`** restricts scoring to the columns that actually differ
between your designs. Everything else holds the same residue in every sequence
and contributes an identical constant to each per-sequence mean, padding the
average and compressing the spread. In the cutinase run, 235 of 276 columns
were identical across all 5,150 designs — 41 varied, of which 38 survived the
gap filter.

Both filters together, on real data:

| Scoring scheme | sd | range / mean |
|---|---|---|
| legacy `mean(A,B)`, all 197 columns | 0.0026 | 4.2% |
| weighted, 38 design columns | 0.0103 | **30.5%** |

Spearman correlation between the two rankings is 0.756, and only **45 of the
legacy top 100 designs remain in the new top 100**. The reordering is real, not
cosmetic.

`--positions` also takes explicit alignment positions (`80,82,94-99`), a file
with one per line, or **a Rosetta resfile** — the last of these needs
`--template` so the residue numbers can be mapped onto columns (see below).
Bare numbers are alignment columns, not residue numbers.

## Calibration

A penalty of 0.2023 is uninterpretable on its own. `--loo` builds a reference
distribution from sequences already known to be real members of the family.

For each natural sequence, a profile is built from **the other** sequences and
the held-out one is scored against it. Scoring the natural set against a
profile built from all of it — what this pipeline originally did — leaks: every
sequence contributed to the distribution it is then measured against, so the
scores come out flatteringly low and cannot serve as a reference.

Each design then gets:

- `percentile_vs_natural` — share of natural folds scoring below it. 0 beats
  every homolog; 100 is worse than all of them. With ~20 folds the resolution
  is about 5 points, so read it as a band.
- `z_vs_natural` — finer ordering, in natural standard deviations.
- `worse_than_worst_natural` — the blunt flag: this design places a residue
  combination no real family member does.

`--loo` also runs a **shuffle control**: each sequence's own residues are
permuted among the scored columns and rescored. This preserves composition
exactly and destroys only residue-to-position pairing, so it isolates whether
the score responds to *where* residues sit. If real sequences don't separate
from their own shuffles, the metric is not measuring position-specific fit and
nothing downstream of it means anything.

On the bundled 22-sequence alignment the control passes comfortably: natural
mean **0.169 ± 0.035** against a shuffle-null mean of **1.123**, a separation of
**27 natural standard deviations**.

A percentile is only meaningful if the reference was built on the same columns
the designs were scored on. Pass matching `--positions` to both stages, or use
`both --loo`, which keeps them consistent by construction — a mismatch is
detected and the annotation skipped rather than silently reported.

## Generating a design space

`make-resfile` runs the logic in reverse. Instead of scoring designs after the
fact, it turns the entropy profile into a Rosetta resfile that constrains what
design is allowed to produce in the first place.

```bash
python entropy_pipeline.py make-resfile \
    --fasta data/natural_cutinase_msa_22.fasta \
    --template "P00590" --template-start 1 --chain A \
    --lock "120,175,188" \
    --compare seq1.resfile \
    --out designed
```

Two decisions, taken from two different parts of the profile:

- **Entropy picks the mode.** `norm_H` at or below `--conserved-max-h` (default
  0.15) means the position is held at the native residue; above it, the
  position is designed.
- **The distribution picks the alphabet.** Amino acids are admitted in
  descending natural probability until the cumulative mass reaches `--top-p`
  (default 0.90). Set sizes then follow each column's own tolerance instead of
  a fixed threshold giving one residue here and fifteen there. The template's
  native residue is always included.

Cysteines are held fixed by default (`--no-lock-cys` to override), since
mutating one half of a disulfide is rarely intended, and `--lock` pins
mechanistically essential residues — a catalytic triad — regardless of entropy.

Output is the resfile plus `*_rationale.csv`, one row per template residue
recording the mode, the alphabet, the entropy and *why* that decision was made.

### Numbering is the thing that will break

A wrong offset fails silently: Rosetta constrains the wrong residue at every
position and reports nothing. So the map is always built from a named sequence
inside the alignment, never assumed — `--template` takes a header substring and
`--template-start` the residue number of its first non-gap residue.

The cutinase case shows why this needs care. The template is in the alignment as
`sp|P00590|CUTI1_FUSVN`, the 230-residue precursor, and matches the 1CEX
deposited sequence at **99.5% identity with a +16 offset** — 1CEX begins at
precursor residue 17, after the signal peptide. So a residue number quoted in
one convention is 16 off in the other, and a catalytic position quoted from a
paper in mature numbering is off by about 31. Confirm `--template-start`
against the ATOM records of the PDB you will actually pass Rosetta; an RCSB
FASTA gives SEQRES order, which need not match author numbering.

### Validate against a resfile you trust

`--compare` diffs the generated specification against an existing one. On the
22-sequence cutinase alignment against a hand-built `seq1.resfile`:

| | |
|---|---|
| reference positions mapped | **193 / 193** |
| same command | 157 / 193 (81%) |
| manual alphabet ⊆ generated | **144 / 157 (92%)** |
| mean alphabet size | manual 2.1, generated 4.3 |

Every reference position mapped, which is the check that the numbering is
right — a bad offset leaves orphans. The generator is roughly twice as
permissive as the manual file at `--top-p 0.90`; lower it to tighten.

The comparison also flags a resfile whose header has a **commented-out default
command** (`#NATRO`). That leaves unlisted residues falling through to Rosetta's
own default, where they are designable to any amino acid — usually the opposite
of what was meant.

### The circularity trap

If you generate a resfile from profile P and then score the resulting designs
with profile P, **the designs score well by construction** — you constrained
them to. The scoring stage stops being a validation and becomes a tautology
that looks like success.

Constraint and evaluation have to draw on different information. Either split
the homologs into disjoint sets, building the resfile from one and calibrating
against the other, or evaluate with something orthogonal: an ESM or ProteinMPNN
likelihood, Rosetta energy, or measured data. Decide this before generating
anything, because it is hard to retrofit.

## Install

```bash
git clone <your-repo-url>
cd entropy-pipeline
```

The core pipeline runs on the **Python standard library alone** (3.9+) — no
install step. `matplotlib` and `numpy` are imported lazily, only when you pass
`--plot`:

```bash
pip install -r requirements.txt   # only for figures / the notebook
```

[MAFFT](https://mafft.cbrc.jp/alignment/software/) is required to produce the
alignments but is not called by this script — you align first, then feed the
result in.

## Quickstart

A 22-sequence natural cutinase alignment ships in `data/`, so this runs with no
setup:

```bash
python entropy_pipeline.py natural \
    --fasta data/natural_cutinase_msa_22.fasta \
    --out natural --loo --plot
```

That writes the reusable profile, the per-position CSV, a calibration file and
a PNG, and prints the shuffle control. Everything in `examples/` was produced
by exactly this command.

## Full workflow

```bash
# 1 — align the natural homologs
mafft --auto natural.fasta > natural_aligned.fasta

# 2 — add designed sequences into that alignment
#     (--add keeps the natural columns fixed, which is what makes
#      column indices comparable between the two stages)
mafft --add designed.fasta --reorder natural_aligned.fasta > combined.fasta

# 3 — both stages, calibrated, scoring only the designable columns
python entropy_pipeline.py both \
    --natural natural_aligned.fasta \
    --combined combined.fasta \
    --design-prefix "design_" \
    --positions auto --loo \
    --out results --plot
```

Then rank on `mean_penalty_weighted` in `results_designed_summary.csv`, and read
`percentile_vs_natural` to see where each design sits against real homologs.

`--design-prefix` is matched case-insensitively against the full FASTA header.
Check the split report before trusting anything — a prefix matching nothing
silently yields zero designs.

## Outputs

| File | Contents |
|---|---|
| `*_profile.json` | Per-column BLOSUM distribution, entropy, `scorable` flag |
| `*_profile.csv` | Same, flat: `position, gap_frac, n_seqs, scorable, H_bits, norm_H, freq_*, blosum_*` |
| `*_summary.csv` | One row per sequence: mean/sum/max/min of each metric, plus calibration fields |
| `*_calibration.json` | Leave-one-out distribution, per-fold scores, shuffle null |
| `*.resfile` | Rosetta design specification (`make-resfile`) |
| `*_rationale.csv` | Per-residue mode, alphabet, entropy and reason |
| `*_scores.json` | Full per-residue scores — **`--per-position` only** |
| `*_per_position.csv` | Long-format per-residue table — **`--per-position` only** |

The last two scale as sequences × columns. The original 5,150-design run wrote a
**319 MB JSON** and an **87 MB CSV**, past what GitHub accepts, so they are now
opt-in. The summary CSV is always written and is what you rank on.

## Notebook

`entropy_pipeline.ipynb` is a cell-by-cell version of the same code, plus a
`logomaker` sequence-logo stage the CLI lacks. Outputs are cleared — run top to
bottom. The config cell near the top sets input filenames and design prefix in
one place. Note that the notebook predates the weighted metric and the column
filters; the CLI is the current implementation.

## Repository layout

```
entropy_pipeline.py        # CLI: natural / designed / both
entropy_pipeline.ipynb     # notebook twin + sequence logos
data/
  natural_cutinase_msa_22.fasta   # aligned, 22 homologs — demo input
  natural_cutinase_23.fasta       # 23-sequence set
  fusarium_msa_27.fasta           # 27-sequence Fusarium alignment
  1cex_reference.fasta            # PDB 1CEX reference sequence
examples/
  profiles/                # Stage 1 + calibration on the 22-seq MSA
  figures/                 # entropy profile, sequence logos
```

## Known limitations

- **The natural set is small.** 22 sequences cluster into 13 groups at 80%
  identity and 9 at 62%, and one pair is 99.6% identical — so roughly 13
  effective observations per column for estimating a 20-way distribution.
  Conservation may reflect shared Fusarium ancestry rather than functional
  constraint. Deepening the alignment (jackhmmer or HHblits against UniRef) is
  the single largest available improvement.
- **No sequence weighting.** That near-duplicate pair votes twice. Henikoff
  position-based weighting would fix it.
- **The softmax is unscaled.** `exp()` over raw BLOSUM scores is very sharp —
  `exp(11)` against `exp(−4)` is a ratio of three million. Standard
  target-frequency pseudocounts (observed counts mixed against BLOSUM62 target
  frequencies, weighted by `N_eff`) would be better behaved and comparable to
  the PSSM literature.
- **Per-column independence.** Coupled positions — salt bridges, packing pairs,
  a catalytic triad acting together — are invisible. A design can score well
  here and still be broken. Mutual information or a DCA model would address it,
  but both need far more effective sequences than 13.
- **No structural term.** Burial, secondary structure and solvent accessibility
  all modulate substitution cost, and none are used.
- **Gaps are skipped, not modeled.** Indel-heavy designs are scored on fewer
  positions than their length suggests; `n_positions` in the summary tells you
  how many.

## Status

Exploratory research code from a working design campaign, cleaned up and made
reproducible. It is not a validated predictor — treat the scores as a triage
filter for prioritizing candidates, not as a prediction of stability or
activity. Nothing here has been checked against measured data, and the method
it implements is established elsewhere (see Relation to existing methods).

A worthwhile exercise for anyone extending this: run PROSS or FuncLib on the
same template and diff their design space against `make-resfile`'s, the same
way `--compare` diffs against a hand-built resfile. Agreement would be real
validation; disagreement would localise exactly what the missing structural
term buys.

## References

### Methods this reimplements

- Goldenzweig A, *et al.* Automated structure- and sequence-based design of
  proteins for high bacterial expression and stability. *Mol Cell* 63, 1–10
  (2016). doi:10.1016/j.molcel.2016.06.012 — **PROSS**
- Weinstein JJ, Goldenzweig A, Hoch SY, Fleishman SJ. PROSS 2: a new server for
  the design of stable and highly expressed protein variants. *Bioinformatics*
  37, 123–125 (2021). doi:10.1093/bioinformatics/btaa1071
- Peleg Y, *et al.* Community-wide experimental evaluation of the PROSS
  stability-design method. *J Mol Biol* 433, 166964 (2021).
  doi:10.1016/j.jmb.2021.166964
- Khersonsky O, *et al.* Automated design of efficient and functionally diverse
  enzyme repertoires. *Mol Cell* (2018). PMID:30270109 — **FuncLib**, the
  closest relative to `make-resfile`
- Weinstein JY, *et al.* Designed active-site library reveals thousands of
  functional GFP variants. *Nat Commun* 14, 2890 (2023).
  doi:10.1038/s41467-023-38099-z — **htFuncLib**
- Goldenzweig A, Fleishman SJ. Principles of protein stability and their
  application in computational design. *Annu Rev Biochem* 87, 105–129 (2018).
  doi:10.1146/annurev-biochem-062917-012102

### Entropy-based design

- Chan CH, *et al.* Relationship between local structural entropy and protein
  thermostability. *Proteins* (2004)
- Bae E, *et al.* Bioinformatic method for protein thermal stabilization by
  structural entropy optimization. (2008)
- Borgo B, Havranek JJ. Automated selection of stabilizing mutations in
  designed and natural proteins. (2012)

### Context

- Isakova *et al.* (2026) — Descent from a common ancestor restricts
  exploration of protein sequence space
- Further reading in the project notes on relative-entropy protein design,
  evolution-entropy hybrid residue ranking, local entropy in proteins, and
  thermal stability in plastic-degrading enzymes

Reference PDFs are not included in this repository (copyrighted). The
entropy-based entries are given by author, title and year only; fill in journal
volume and pages from the PDFs before reusing these in a manuscript.

## License

MIT — see `LICENSE`.
