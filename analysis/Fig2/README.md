## Script

`fig2.m` — run as a whole or section by section.

## Inputs

| File | Contents |
| --- | --- |
| `data/mouse/Scores/Virulence_screen_clean_table.csv` | Mono-colonization screen, one row per mouse per day |
| `data/mouse/Scores/ProtectionScreen_CDI_mouse.csv` | Co-colonization (ST1 + VPI) screen, same schema |

## What it computes

The two tables are concatenated and a per-mouse key `exp_id = experiment_strain_mouse` is built, then split into two arms:

- **Virulence (mono-colonization).** `relweight` is sign-flipped so that a larger score means more weight loss, i.e. more virulence. Model: `relweight ~ cdiffstrain + (1|day) + (1|exp_id)` with `fitlme`, reference level `ui`.
- **Protection (co-colonization).** reference level `vpi`. Here a larger fixed effect means less weight loss than VPI alone, i.e. more protection.

Fixed effects and 95% CIs from `coefCI` become the per-strain virulence and protection scores. The relationship between the two scores is tested with Spearman (and Pearson) correlation and a `fitlm` trend line.

Survival curves use `ecdf(..., 'Censoring', ~death, 'Function','survivor')` on the last observed day per mouse.

## Outputs

Four figure windows:

1. Weight trajectories for VPI, ST1-75 + VPI and ST1-68 + VPI (mean ± SD across experiments).
2. Survival curves for the same three arms.
3. Protection score per strain, sorted, with 95% CIs.
4. Mono-colonization virulence vs. co-colonization protection, ST1-75 and ST1-68 highlighted, with Spearman ρ and p.

