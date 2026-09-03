# Figure 3 — BIOLOG substrate use, niche breadth and predicted outcomes

Two things live here: the **preprocessing pipeline** that turns raw BIOLOG PM1 kinetic plates into binary growth calls, and the **figure script** that relates those calls to the protection scores from Figure 2.

## Scripts

| Script | Role |
| --- | --- |
| `Fig3.m` | Figure panels. Reads the committed growth matrix. |
| `run_biolog_folders.m` | Preprocessing driver: loops over strain folders, calls the plate analyzer, writes per-strain growth calls. |
| `analyze_biolog_plate.m` | Scores one PM1 plate: mixed model per well, Bonferroni correction, QC on controls. |
| `combine_binary_growth.m` | Merges per-strain calls into `data/Biolog/Biolog_growth_matrix.xlsx`. |
| `biologDataDir.m` | Helper returning the absolute path to `data/Biolog`. |

## Inputs

| File | Contents |
| --- | --- |
| `data/Biolog/Biolog_growth_matrix.xlsx` | 96 substrates × 21 ST1 strains + `VPI`, binary (`1` = growth) |
| `data/Biolog/strain_groups.xlsx` | Per-strain `Protection_Estimate`, `Virulence_Estimate`, groups, `Total_Metabolites` |
| `data/Biolog/moas.xlsx` | BIOLOG chemical class mapping |
| `data/Biolog/Biolog_PM1_files/` | Raw plates and per-strain calls |

##

How a plate is scored (`analyze_biolog_plate.m`):

1. Read the plate (`Duration_Hours_` plus 96 wells), reshape to long form, keep timepoints `< 24 h` (`TimeCutoff`).
2. Join well IDs to substrate names via `Biolog_names_sorted.xlsx`.
3. Baseline-subtract each well at its own first timepoint, flooring at 0.
4. Fit `AbsorbanceminusT0 ~ Time*Names + (Time|Well)` — growth is the substrate-specific slope relative to `Negative Control`.
5. Call growth where the interaction term is positive and `p ≤ Alpha/NumTests` (Bonferroni, default `0.05/95`).
6. QC: `a-D-Glucose` must be called `1`, `Negative Control` must be called `0`; failures raise warnings, they do not stop the run.

Strains with replicate plates (`*_PM1_1.xlsx`, `_2`, `_3`) get one call per replicate plus a consensus call that is the **logical OR across replicates** — a substrate counts as used if any replicate supports growth. `combine_binary_growth` ingests only strain-level `*_binary_growth.xlsx` files, skipping the numbered replicate files and the combined output itself.


## `Fig3.m`

Strains are sorted by protection score; substrates are ordered by VPI use, then chemical class, then prevalence across strains.

- **Substrate-use matrix** — strains (sorted by protection) × substrates, with a vertical rule separating substrates VPI does and does not use.
- **Pairwise resource classes** — shared / ST1-private / VPI-private substrate counts per strain, stacked.
- **Breadth vs. protection** — total substrates used per strain vs. protection score, Spearman ρ and least-squares line.
- **Model outcome classes** — the Figure 1 invasion rule applied to every strain over `λ ∈ [0, 82]` (501 points); bars give the fraction of the λ grid falling in each of the four outcome classes.
- **Chemical-class breadth vs. protection** — per class, Spearman correlation between the number of substrates used in that class and protection, with Benjamini–Hochberg FDR (`*` p < 0.05, `**` q < 0.05).

## Outputs

Five figure windows. Console output: the breadth-vs-protection correlation and a per-class table of ρ, p and q.
