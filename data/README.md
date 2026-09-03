# `data/`

All raw and derived data used by the scripts in [`../analysis/`](../analysis). Everything here is read-only for the figure scripts; the one exception is the BIOLOG preprocessing described below.

```
data/
├── Biolog/            BIOLOG PM1 carbon-source plates and binary growth calls
├── mouse/             Mouse weight/survival screens and RAG1 KO rechallenge
├── qPCR/              1:5 mixed-infection weights and strain fractions
├── 16s_sequencing/    Residual-microbiota genus abundances
├── Genomics/          Gene presence/absence matrix and phylogenies
└── flow_cytometry/    Exported flow events, ungated and per gate
```

**Strains.** 21 ST1 clinical isolates (`ST1-2`, `ST1-6`, `ST1-11`, `ST1-12`, `ST1-19`, `ST1-20`, `ST1-23`, `ST1-25`, `ST1-27`, `ST1-49`, `ST1-53`, `ST1-57`, `ST1-58`, `ST1-62`, `ST1-63`, `ST1-65`, `ST1-66`, `ST1-67`, `ST1-68`, `ST1-69`, `ST1-75`) plus the toxigenic reference `VPI10463` (`VPI`). 
