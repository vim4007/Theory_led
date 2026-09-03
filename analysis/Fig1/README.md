## Script

`Fig1.m` — run as a whole or section by section.

## Inputs

| File | Used for |
| --- | --- |
| `data/Biolog/Biolog_growth_matrix.xlsx` | Binary substrate use per strain; columns `ST1_75`, `ST1_68`, `VPI` |
| `data/Biolog/strain_groups.xlsx` | Loaded for convenience; not used by the current panels |

## What it computes

For a strain with binary use vector `s` and VPI use vector `vpi`:

- `n_shared = sum(s & vpi)` — substrates both strains use
- `n_ST1 = sum(s & ~vpi)` — ST1-private substrates
- `n_VPI = sum(~s & vpi)` — VPI-private substrates

Invasion scores as a function of the resource-limitation threshold λ:

```
I_ST1(λ) = n_ST1 − λ · n_VPI / (n_shared + n_VPI)
I_VPI(λ) = n_VPI − λ · n_ST1 / (n_shared + n_ST1)
```

The sign pair `(I_ST1, I_VPI)` classifies the outcome: coexistence, ST1 excludes VPI, VPI excludes ST1, or neither invades.

## Outputs

Three figure windows, plus the shared/private substrate counts for ST1-75 and ST1-68 printed to the console.

1. Outcome rule — the four quadrants of the `(I_ST1, I_VPI)` plane.
2. Invasion scores vs. λ for ST1-75 against VPI10463-like counts.
3. The same for ST1-68.
