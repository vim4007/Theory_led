## Script

`fig4.m` — run as a whole or section by section.

## Inputs

| File | Contents |
| --- | --- |
| `data/Biolog/Biolog_growth_matrix.xlsx` | Binary substrate use; `Negative Control` dropped |
| `data/qPCR/weights.xlsx` | Relative weights, days 0–3: columns `ST1_75`, `VPI`, and three `mix` replicates |
| `data/qPCR/qPCR section.xlsx` | Per-mouse strain fractions (%) by day: `Mouse`, `Day`, `ST1-75`, `VPI` |

## What it computes

**Predicted invasion robustness.** Shared / ST1-private / VPI-private substrate counts for ST1-75 and ST1-68, then the ST1 invasion score

```
I_ST1(λ) = n_ST1 − λ · n_VPI / (n_shared + n_VPI)
```

**ODE simulation.** A five-variable consumer–resource system (two strain densities, three resource pools) integrated with `ode45` over `t ∈ [0, 25]`:

```
dN_i/dt = N_i (e·u·R_shared + e·u·R_private,i − d − D)
dR_k/dt = n_k·δ0 − D·R_k − u·R_k·(consumers of k)
```

Parameters: `e = 0.5` (conversion efficiency), `d = 0.05` (death), `D = 0.2` (dilution), `δ0 = 0.08` (per-substrate supply), `u = 0.25` (uptake). Substrate counts `n_shared`, `n_ST1`, `n_VPI` are ST1-75's measured values; resources start at their no-consumer steady state `n·δ0/D`, and the strains start at `N_ST1 = 1`, `N_VPI = 5` to mirror the 1:5 inoculum. 

**Experiment.** Relative weight over days 0–3 for ST1-75 alone, VPI alone, and the mix (mean ± SD of three replicates); and mean qPCR strain fractions at days 0.5, 1 and 3 for co-infected mice (mice 3, 4, 5), stacked VPI over ST1-75.

## Outputs

Four figure windows:

1. Predicted invasion score vs. λ, ST1-75 vs. ST1-68.
2. Simulated ST1-75 and VPI densities from a 1:5 start.
3. Measured relative weight, three arms.
4. Measured qPCR composition over time.
