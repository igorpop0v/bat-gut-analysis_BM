# PICRUSt2 Functional Prediction Analysis

## Scope

PICRUSt2 was used to infer the functional potential of the bacterial communities from 16S rRNA gene amplicon data. The resulting gene-family and MetaCyc pathway profiles represent computational predictions based on phylogenetic placement and reference genomes.

Throughout the project, these results are described as **predicted**, **inferred**, or **putative functional potential**. They must not be interpreted as direct measurements of gene expression, metabolic activity, metabolite production, or pathway flux.

## Experimental groups

The analysis included 76 independent biological samples divided into three groups:

- **ABH**: active before hibernation;
- **H**: hibernation;
- **AAH**: active after hibernation.

## Input data

PICRUSt2 received:

1. Deblur representative ASV sequences trimmed to 190 nucleotides;
2. the corresponding unrarefied ASV count table containing biological samples only.

The unrarefied table was used because PICRUSt2 requires the original observed sequence abundances to estimate community functional profiles. Rarefaction was not applied before functional prediction or differential association analysis.

Input preparation is implemented in:

- `workflows/picrust2/scripts/01_prepare_inputs.sh`

## Software environment

PICRUSt2 was run in Docker using the following fixed container image:

```text
quay.io/biocontainers/picrust2:2.6.3--pyhdfd78af_2
```

The pipeline is implemented in:

- `workflows/picrust2/scripts/02_run_picrust2.sh`

The main command uses `picrust2_pipeline.py` with the representative sequences, ASV count table, and a user-specified number of processor threads.

The standard PICRUSt2 workflow includes:

1. phylogenetic placement of study ASVs;
2. hidden-state prediction of gene-family copy numbers;
3. normalization for predicted 16S rRNA gene copy number;
4. prediction of community EC abundances;
5. inference of MetaCyc pathway abundances.

The default PICRUSt2 pathway-inference settings, including MinPath and pathway gap filling, were retained.

## PICRUSt2 quality control

Quality control is implemented in:

- `workflows/picrust2/r/01_picrust2_qc.R`

The following results were obtained:

| Metric | Result |
|---|---:|
| Input ASVs | 446 |
| Successfully placed ASVs | 445 |
| Unplaced ASVs | 1 |
| Successfully placed ASVs | 99.78% |
| Bacterial placements | 445 |
| Archaeal placements | 0 |
| Mean ASV NSTI | 0.0676 |
| Median ASV NSTI | 0.0544 |
| Maximum ASV NSTI | 1.0091 |
| ASVs with NSTI > 0.15 | 28 |
| ASVs with NSTI > 0.30 | 4 |
| ASVs with NSTI > 1 | 1 |
| ASVs with NSTI > 2 | 0 |
| Mean sample-weighted NSTI | 0.0617 |
| Median sample-weighted NSTI | 0.0614 |
| Minimum sample-weighted NSTI | 0.0161 |
| Maximum sample-weighted NSTI | 0.1032 |

The sample-weighted NSTI distributions did not differ significantly among the experimental groups:

```text
Kruskal-Wallis p = 0.712
```

Consequently, there was no evidence that prediction quality differed systematically among ABH, H, and AAH samples.

NSTI was retained as a quality-control measure and was not included as a covariate in the primary biological models. Weighted NSTI is derived from microbial composition, and including it in the model could partially adjust away the biological signal of interest.

## Preparation of R input files

PICRUSt2 results were converted into analysis-ready R input files by:

- `workflows/picrust2/scripts/03_prepare_r_inputs.sh`

The primary functional table contains unstratified predicted MetaCyc pathway abundances. Unstratified profiles were used because the primary objective was to compare community-level predicted functional potential among the three experimental groups.

## Predicted pathway composition

Predicted pathway composition is analysed in:

- `workflows/picrust2/r/02_pathway_composition.R`

Predicted pathway abundances were converted to relative abundances by total-sum scaling separately within each sample:

```text
relative abundance = pathway abundance / total pathway abundance in the sample
```

Mean relative abundance was then calculated within ABH, H, and AAH.

The primary composition heatmap contains 25 pathways. Pathways were ranked first by their maximum mean relative abundance across the three groups and then by their overall mean relative abundance across all samples.

The top-25 display is a visualization choice only. All predicted pathways remain available in the complete exported tables.

## Predicted pathway diversity

Predicted pathway diversity is analysed in:

- `workflows/picrust2/r/03_pathway_diversity.R`

### Normalization

The complete predicted pathway table was converted to relative abundance using total-sum scaling. Rarefaction was not used because the pathway values are inferred functional abundances rather than directly observed sequencing counts.

### Alpha diversity

The primary alpha-diversity metrics were:

- Shannon index;
- Simpson index;
- Pielou evenness.

Overall differences among groups were tested using the Kruskal-Wallis test. Pairwise comparisons were performed using Wilcoxon rank-sum tests, with Benjamini-Hochberg correction applied independently within each diversity metric.

The principal results were:

| Metric | Overall result | Significant pairwise result |
|---|---|---|
| Shannon | p = 0.0245 | ABH vs H, q = 0.0249 |
| Simpson | p = 0.251 | none |
| Pielou evenness | p = 0.0050 | H vs AAH, q = 0.0044 |

These results indicate changes in predicted functional diversity and evenness during hibernation, although the specific conclusion depends on the diversity metric.

### Beta diversity

Bray-Curtis dissimilarity was calculated from TSS-normalized predicted pathway profiles.

Overall group differences were tested using PERMANOVA with 9,999 permutations:

```text
pseudo-F = 11.15
R² = 0.234
p < 0.001
```

Pairwise PERMANOVA results were:

| Comparison | q-value |
|---|---:|
| ABH vs H | 0.00015 |
| ABH vs AAH | 0.349 |
| H vs AAH | 0.00015 |

PERMDISP was used to evaluate the homogeneity of multivariate dispersion:

```text
F = 2.76
p = 0.066
```

No pairwise PERMDISP comparison remained significant after Benjamini-Hochberg correction.

Therefore, the significant PERMANOVA result is most consistent with differences in group centroids rather than a clear difference in within-group dispersion.

The beta-diversity results indicate that the predicted pathway profile during hibernation differed from both active states, whereas ABH and AAH were not significantly different at the community-wide pathway level.

## Differential abundance and prevalence analysis

Differential association analysis is implemented in:

- `workflows/picrust2/r/04_maaslin3_pathways.R`

MaAsLin3 was selected because it separately models:

1. **abundance**: how much of a pathway is predicted among samples where it is detected;
2. **prevalence**: the probability that a pathway is detected in a sample.

### Filtering

A pathway was retained if it was detected in at least 10% of samples.

With 76 samples, this required detection in at least eight samples.

| Filtering stage | Pathways |
|---|---:|
| Before filtering | 470 |
| After filtering | 439 |
| Removed | 31 |

The filter removed very rare pathways with insufficient information for stable group comparisons while retaining most predicted functions.

### Model

The primary MaAsLin3 model was:

```text
~ group + log2_library_size
```

The original unrarefied 16S sequencing depth was included as `log2_library_size`. This is particularly important for prevalence models because deeper sequencing can increase the probability of detecting low-abundance features.

The analysis used:

- H as the reference group;
- TSS normalization;
- base-2 log transformation;
- augmented logistic regression for prevalence;
- linear regression for non-zero abundance;
- standardization of continuous covariates;
- median coefficient comparison for abundance;
- a zero null hypothesis for prevalence;
- Benjamini-Hochberg false-discovery-rate correction;
- an individual significance threshold of `q < 0.05`.

The primary comparisons were:

- ABH vs H;
- AAH vs H.

Positive coefficients indicate a higher predicted pathway abundance or prevalence in the corresponding active group relative to H. Negative coefficients indicate a lower value in the active group and, consequently, a higher value during hibernation.

### MaAsLin3 results

The primary analysis identified:

| Result | Number |
|---|---:|
| Significant association rows | 384 |
| Unique pathways with at least one significant association | 210 |
| Pathways significant in both group comparisons | 147 |

The significant results were distributed as follows:

| Comparison | Abundance | Prevalence |
|---|---:|---:|
| ABH vs H | 131 | 52 |
| AAH vs H | 144 | 57 |

A single pathway can produce more than one significant result because it can be associated with both group comparisons and with both abundance and prevalence.

The large number of significant results does not imply an equivalent number of independent biological mechanisms. MetaCyc pathways are highly correlated, share reactions and enzymes, and are inferred from the same underlying taxonomic profiles.

## Custom MaAsLin3 visualization

The custom visualization is implemented in:

- `workflows/picrust2/r/05_visualize_maaslin3_pathways.R`

The primary figure displays 25 unique pathways.

Pathways were ranked by:

1. the smallest individual q-value across both group comparisons and both models;
2. the maximum absolute difference between the coefficient and its null hypothesis, used only to resolve ties.

For every selected pathway, all available ABH vs H and AAH vs H abundance and prevalence results are displayed.

The figure uses:

- circles for abundance associations;
- triangles for prevalence associations;
- horizontal lines for 95% confidence intervals;
- purple shading for abundance q-values;
- green shading for prevalence q-values;
- white symbols for `q ≥ 0.05`;
- separate null-hypothesis lines for abundance and prevalence.

Values below `0.000001` are displayed as `≤0.000001` to avoid infinite legend labels when numerical p-values or q-values underflow to zero.

The sequencing-depth covariate is retained in the statistical model but is not displayed in the primary biological figure.

The complete results remain available in TSV format. The top-25 figure is not a replacement for the full result table.

## Biological interpretation

The predicted functional profiles indicate substantial remodelling of bacterial community functional potential during hibernation.

Relative to H, both active groups showed higher predicted representation of several pathways associated with:

- peptidoglycan and cell-envelope biosynthesis;
- O-antigen building-block biosynthesis;
- carbohydrate utilization;
- degradation of several organic compounds.

Several pathways associated with anaerobic metabolism, fermentation, nucleotide salvage, amino-acid degradation, and cofactor biosynthesis showed the opposite direction and were predicted to be higher during hibernation.

The broadly similar directions of ABH vs H and AAH vs H associations are consistent with the beta-diversity result in which H differed from both active groups, whereas ABH and AAH were not significantly different.

However, the current pathway-level MaAsLin3 model did not directly test ABH vs AAH. Similar effect directions must therefore not be interpreted as formal equivalence between these groups.

The functional results are expected to reflect, at least partly, the taxonomic changes observed in the 16S analysis because PICRUSt2 predictions are calculated from the abundance and phylogenetic placement of the same ASVs. The functional analysis therefore provides a biologically useful interpretation of the taxonomic shifts but is not independent confirmation of them.

## Important limitations

### Predicted potential is not measured activity

PICRUSt2 predicts potential gene and pathway content. It does not measure:

- gene expression;
- RNA abundance;
- protein abundance;
- metabolite concentration;
- pathway flux;
- actual physiological activity.

Terms such as “increased pathway activity” must not be used when describing these results.

### Reference-genome dependence

Predictions depend on the availability and quality of related reference genomes. Strain-level gene-content differences and horizontal gene transfer cannot be reliably inferred from a short 16S marker.

### Pathway labels may be overly specific

Pathway names may be inferred from shared enzymes and partial pathway evidence. Consequently:

- `methanogenesis from acetate` is not evidence that methanogens were present, particularly because no archaeal placements were observed;
- a pathway labelled `β-lactam resistance` is not phenotypic evidence of antimicrobial resistance;
- predicted toluene or catechol degradation pathways are not evidence of exposure to those compounds.

Such results may be reported as predicted pathway signals but should not be interpreted literally without independent genomic, transcriptomic, metabolomic, or culture-based validation.

### Correlated pathways

Many MetaCyc pathways share reactions and predicted gene families. Significant pathways are therefore not statistically or biologically independent.



## Reproducibility

The complete workflow is represented by the following scripts:

```text
workflows/picrust2/scripts/01_prepare_inputs.sh
workflows/picrust2/scripts/02_run_picrust2.sh
workflows/picrust2/scripts/03_prepare_r_inputs.sh

workflows/picrust2/r/01_picrust2_qc.R
workflows/picrust2/r/02_pathway_composition.R
workflows/picrust2/r/03_pathway_diversity.R
workflows/picrust2/r/04_maaslin3_pathways.R
workflows/picrust2/r/05_visualize_maaslin3_pathways.R
```

Large input files, Docker outputs, and generated analysis results are stored outside the Git repository. The repository contains code, configuration templates, methodological decisions, and documentation required to reproduce the analysis.

R package and system versions are recorded in the generated `session_info.txt` files.

## References

1. Douglas GM, Maffei VJ, Zaneveld JR, et al. PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology*. 2020;38:685–688. [https://doi.org/10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6)

2. Nickols WA, Kuntz T, Shen J, et al. MaAsLin 3: refining and extending generalized multivariable linear models for meta-omic association discovery. *Nature Methods*. 2026;23:554–564. [https://doi.org/10.1038/s41592-025-02923-9](https://doi.org/10.1038/s41592-025-02923-9)

3. PICRUSt2 documentation: Key Limitations. [https://github.com/picrust/picrust2/wiki/Key-Limitations](https://github.com/picrust/picrust2/wiki/Key-Limitations)

4. MaAsLin3 Bioconductor tutorial. [https://bioconductor.posit.co/packages/3.21/bioc/vignettes/maaslin3/inst/doc/maaslin3_tutorial.html](https://bioconductor.posit.co/packages/3.21/bioc/vignettes/maaslin3/inst/doc/maaslin3_tutorial.html)
