# 16S rRNA analysis in R

The scripts are intended to be run in numerical order:

1. `01_parse_taxonomy.R` — parse and clean GTDB taxonomy;
2. `02_microbiota_composition.R` — summarize and visualize community composition;
3. `03_diversity_analysis.R` — analyze alpha and beta diversity;
4. `04_maaslin3_daa.R` — differential abundance analysis with MaAsLin3.

Input and output directories are defined in `local_paths.R`.
Create this local file from `local_paths.R.example`.

The local paths file, input data, and generated results are not tracked by Git.
