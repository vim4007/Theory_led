## Scripts

`Fig5.m` is the script to run: it builds every panel of the figure and carries the shared helper functions (`flow_fractions`, `plot_immune_bars`, `bh_fdr`).


## Inputs

| File / folder | Panel | Contents |
| --- | --- | --- |
| `data/flow_cytometry/ks10_adaptive_csv_files/` | A | Per-mouse event tables: `all_adaptive/` |
| `data/flow_cytometry/ks10_innate_csv_files/` | B | Same layout: `all_cells/` plus one folder per gate |
| `data/mouse/rag1ko/KS11_rechallenge_relweight.xlsx` | C | `day`, `mouse_id`, `group`, `weight_g`, `relweight` |
| `data/16s_sequencing/tblAbund.xls` | D | 
| `data/Genomics/binary_gene_pre_abs.csv` | E | 21 strains × 5,109 genes (binary), plus `Virulence_Estimate` and `Protection_Estimate` |

## What it computes

**A / B — immune fractions.** 

**C — RAG1 KO rechallenge.** Mean ± SEM relative weight per day for `uninfected_b6`, `st175_b6` and `st175_rag1ko`.

**D — residual microbiota.** `Clostridioides` is removed so that only the residual community is scored. Per sample, Shannon diversity of the renormalized residual genera. Restricted to `mnvc` antibiotic treatment, days 0, 1, 2, 7, comparing ST1-75 vs. ST1-12 initial colonization. Three test families, each BH-corrected separately:

1. Shannon per day (rank-sum);
2. Bray–Curtis group separation per day (5,000-permutation test, `rng(42)`);
3. the top-12 genera by mean abundance at day 1 (rank-sum).

**E — genomics.**  PCA on the mean-centred accessory matrix, points coloured by protection score. `genomic_analysis.m`.

## Outputs

Five figure windows from `Fig5.m`. 
