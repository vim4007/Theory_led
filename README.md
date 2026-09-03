# Metabolic niche breadth determines protective competition among *Clostridioides difficile* strains
## Overview
Antibiotic disruption of colonization resistance can allow toxigenic *Clostridioides difficile* to expand in the gut and cause disease. Hosts colonized with non-toxigenic or weakly virulent *C. difficile* strains can have a lower risk of disease, but which strains protect against a more virulent one remains difficult to predict. Here, we propose that ecological niche breadth provides a quantitative framework for identifying protective strains. We developed a consumer-resource model in which a candidate protective strain and toxigenic VPI10463 compete for shared and strain-private resource classes. The model predicts protection when the candidate strain can persist with or exclude VPI10463 under resource-limiting conditions. In mouse experiments, genetically related sequence type 1 (ST1) strains isolated from MSKCC patients differed widely in their ability to protect against disease when co-colonized with VPI10463. Protection was inversely associated with mono-colonization virulence of ST1 strains, but low virulence alone did not guarantee protection. BIOLOG substrate-use profiles showed that protective strains tended to use more carbon sources overall, including substrates that VPI10463 did not use. An experimental test of protective ST1-75 supported the model: ST1-75 increased from a 1:5 starting disadvantage to near dominance during mixed infection with VPI10463. Immune, residual-microbiome, and genomic analyses did not identify a single alternative strain-level explanation for protection. These data support metabolic niche breadth as a measurable ecological trait that can help rank candidate *C. difficile* strains for protection against toxigenic competitors.


## Repository structure

```
analysis/          Analysis code, one folder per figure, each with its own README
data/              All raw and derived data, documented in data/README.md
figures/           Final assembled figures (.ai, .pdf) and figure legends
Supplementary/     Supplementary figure scripts and assembled panels
```

- [`analysis/README.md`](analysis/README.md) — how the code is organized, conventions, and how to run it
- [`data/README.md`](data/README.md) — what every data file contains and where it came from

## Reproducibility

```matlab
addpath(genpath('/path/to/Theory_led'));
Fig1    % then fig2, Fig3, fig4, Fig5
```

Scripts locate the repository root relative to their own location, so no paths need editing — but the repository must be on the MATLAB path and scripts must stay at `analysis/<FigN>/`. Each script is a cell script: run the whole file, or step through it section by section to build one panel at a time. Panels open as figure windows and summary statistics are printed to the console; final figures were exported and assembled in Illustrator.

Data files are read-only for the figure scripts. The single exception is the BIOLOG preprocessing in `analysis/Fig3/`, which regenerates `data/Biolog/Biolog_growth_matrix.xlsx`; its output is already committed, so it does not need to be re-run.

## Dependencies

MATLAB R2022b or newer
Required toolboxes: Statistics and Machine Learning Toolbox

## Citation

## Contact
Vishwas Mishra, vim4007@med.cornell.edu
