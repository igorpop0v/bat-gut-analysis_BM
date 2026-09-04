# Bat gut microbiome, predicted functional potential, and virome analysis

Reproducible workflows for the analysis of faecal samples collected from common noctules (*Nyctalus noctula*) before, during, and after hibernation in a wildlife rehabilitation setting.

The repository supports three complementary parts of the project:

1. bacterial and archaeal community profiling using 16S rRNA gene amplicon sequencing;
2. prediction of community functional potential with PICRUSt2;
3. faecal RNA virome profiling, including bacteriophage summaries and comparative analysis of densovirus genomes.

Only code, small configuration tables, metadata required to interpret the code, and methodological documentation are version-controlled. Raw reads, reference databases, QIIME 2 artifacts, intermediate files, result tables, and figures are intentionally excluded.

## Study design

The following group abbreviations are used throughout the repository:

- **ABH** — active bats before hibernation, sampled shortly after admission from the wild;
- **H** — bats sampled during hibernation in rehabilitation;
- **AAH** — active bats sampled after hibernation and before release.

The 16S and RNA-virome datasets were generated from different faecal specimens. They are therefore analysed as complementary datasets rather than as paired multi-omics measurements.

## Repository structure

```text
workflows/
├── 16s/        # QIIME 2 and R analyses of 16S rRNA data
├── picrust2/   # PICRUSt2 execution, QC, and R analyses
└── virome/     # RNA-virome, bacteriophage, assembly, and densovirus analyses

docs/
├── 16s/        # 16S QC, diversity, and differential-association decisions
├── picrust2/   # interpretation and analysis of predicted functions
└── virome/     # virome sensitivity and densovirus workflow documentation
```

Scripts are numbered in their intended execution order within each workflow. The virome workflow contains both community-level analyses and a longer genome-level densovirus branch.

## 16S rRNA workflow

The amplicons target the V3–V4 region. Because reverse-read quality was inadequate for reliable paired-end reconstruction, the primary workflow uses primer-trimmed forward reads. Reads are quality-filtered and denoised with Deblur at a fixed length of 190 nucleotides. Taxonomy is assigned with the GTDB r232 classifier, and interpretation is limited to the genus level.

The shell scripts in `workflows/16s/scripts/` perform:

1. FastQC and MultiQC assessment;
2. forward-read manifest generation and QIIME 2 import;
3. primer removal and quality filtering;
4. Deblur denoising and artifact summaries;
5. taxonomic classification and table export;
6. rarefaction, alpha diversity, and Bray–Curtis beta diversity;
7. export of stable input tables for R.

The R scripts in `workflows/16s/r/` parse taxonomy, visualize microbial composition, analyse diversity, and perform MaAsLin3 abundance and prevalence modelling.

The biological-sample diversity analysis uses a rarefaction depth of 10,000 reads. The positive control is retained for technical QC but excluded from biological comparisons.

Key decisions are documented in:

- [`docs/16s/qc_decision.md`](docs/16s/qc_decision.md)
- [`docs/16s/diversity_metrics.md`](docs/16s/diversity_metrics.md)
- [`docs/16s/differential_association_analysis.md`](docs/16s/differential_association_analysis.md)

## PICRUSt2 workflow

PICRUSt2 receives the unrarefied biological-sample ASV table and the matching Deblur representative sequences. The workflow predicts gene-family and MetaCyc pathway profiles, assesses phylogenetic placement and NSTI, summarizes pathway composition, analyses pathway diversity, and performs MaAsLin3 modelling.

These outputs describe **predicted functional potential**. They are not direct measurements of transcription, metabolic activity, metabolites, or pathway flux.

The complete rationale, QC results, normalization choices, statistical strategy, and interpretation limits are described in [`docs/picrust2/function_prediction_analysis.md`](docs/picrust2/function_prediction_analysis.md).

## RNA-virome workflow

The initial virome scripts in `workflows/virome/scripts/` perform:

1. FastQC/MultiQC and paired-read trimming with fastp;
2. removal of reads matching the host reference;
3. removal of ribosomal RNA with SortMeRNA;
4. taxonomic classification of non-rRNA reads with MetaBuli;
5. preparation of viral-family and bacteriophage-family count tables.

The R scripts in `workflows/virome/r/` generate viral-composition figures, primary and sensitivity diversity analyses, MaAsLin3 models, and a clearly labelled exploratory MaAsLin2 method-sensitivity analysis.

The primary virome analysis retains all ten libraries and is not rarefied. This preserves the very small hibernation group but requires explicitly exploratory interpretation. Sensitivity analyses evaluate the effects of the low-depth sample and the dominant *Parvoviridae* signal. Details are provided in [`docs/virome/diversity_sensitivity_analysis.md`](docs/virome/diversity_sensitivity_analysis.md).

## Densovirus genome analysis

The genome-level virome branch includes per-sample MEGAHIT assembly, QUAST summaries, MetaBuli classification of contigs, recovery of Parvoviridae candidates, read mapping, depth assessment, comparison with earlier study genomes and RefSeq NC_031450.1, orientation standardization, protein-marker phylogeny, whole-genome similarity, and conservative feature annotation.

The final phylogenetic analysis uses partitioned NSP1 and SP1 amino-acid alignments. MAFFT alignments are trimmed by column occupancy and analysed with IQ-TREE; primary, study-inclusive, and sensitivity trees are retained to document robustness. The manuscript-facing tree is the study-inclusive tree because it retains all study densovirus sequences with sufficient marker information.

Only read-supported and protein-validated CDS features are written to provisional annotations. Unresolved ORFs are reported rather than forced into GenBank-ready feature calls. The generated submission materials therefore remain provisional until final manual review.

See [`docs/virome/densovirus_genome_analysis.md`](docs/virome/densovirus_genome_analysis.md) for the full strategy and interpretation boundaries.

## Reproducing the analyses

Docker is required for the command-line workflows. R 4.x is used for statistical analysis and plotting.

Clone the repository and enter it:

```bash
git clone https://github.com/igorpop0v/bat-gut-analysis_BM.git
cd bat-gut-analysis_BM
```

Create local path files from the supplied examples and edit only the copies:

```bash
cp workflows/16s/r/local_paths.R.example workflows/16s/r/local_paths.R
cp workflows/picrust2/r/local_paths.R.example workflows/picrust2/r/local_paths.R
cp workflows/virome/r/local_paths.R.example workflows/virome/r/local_paths.R
```

The local path files are ignored by Git. This keeps machine-specific absolute paths out of the public repository.

Core R dependencies include `dplyr`, `tidyr`, `readr`, `tibble`, `stringr`, `purrr`, `ggplot2`, `patchwork`, `vegan`, `qiime2R`, `maaslin3`, `Maaslin2`, `ggnewscale`, `ggrepel`, `ggtext`, `pheatmap`, `gridExtra`, `ape`, and `ggtree`.

Pinned software and container versions are listed in [`docs/software_versions.md`](docs/software_versions.md).

## Data availability

The 16S rRNA dataset is associated with NCBI BioProject [PRJNA1175038](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1175038). RNA-sequencing records and curated densovirus genome submissions will be linked after deposition and accession assignment.

Raw reads and large derived files are not duplicated in this repository. Small public metadata files are included when they are required to reproduce group definitions, filtering, labelling, or reference selection.

## Interpretation note

The RNA-virome dataset is small and highly heterogeneous in viral read recovery. Community-level virome statistics are exploratory. Detection of densovirus sequences in DNase-treated faecal RNA libraries can reflect viral transcripts or protected viral nucleic acid associated with ingested arthropod material; it does not by itself demonstrate productive infection of the bat host.
