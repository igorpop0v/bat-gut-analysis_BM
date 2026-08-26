# Differential abundance and prevalence analysis

## Purpose

Differential association analysis was used to identify bacterial taxa whose
abundance or prevalence differed among the three physiological groups of bats:

- ABH: active before hibernation;
- H: hibernation;
- AAH: active after hibernation.

The analysis was performed with MaAsLin3, which models abundance and
prevalence associations separately.

## Samples and study design

The analysis included 76 biological samples.

Hibernation was used as the reference group. The primary contrasts were:

- ABH versus H;
- AAH versus H.

A positive coefficient therefore indicates a higher abundance or prevalence
in ABH or AAH relative to hibernation.

## Input data

The analysis used the ASV count table produced after Deblur
denoising.

The ASV counts were joined with the parsed GTDB taxonomy and independently
aggregated at the following taxonomic ranks:

- phylum;
- class;
- order;
- family;
- genus.

Species-level analysis was not performed because the 190-nucleotide 16S rRNA
marker fragment does not provide sufficient resolution for definitive
species-level assignments.

## Taxon filtering

Taxa were retained if they were detected in at least 10% of the 76 biological
samples.

This corresponds to a minimum of 8 positive samples:

```text
ceiling(0.10 × 76) = 8
```

The filtering results were:

| Taxonomic rank | Before filtering | After filtering | Removed |
|---|---:|---:|---:|
| Phylum | 16 | 9 | 7 |
| Class | 20 | 13 | 7 |
| Order | 52 | 36 | 16 |
| Family | 86 | 57 | 29 |
| Genus | 137 | 84 | 53 |

The prevalence threshold was applied independently at each taxonomic rank.

The purpose of this filtering was to:

- remove extremely sparse taxa for which model estimates would be unstable;
- reduce the number of poorly supported prevalence models;
- reduce the multiple-testing burden;
- improve the interpretability of the results.

The threshold was selected before evaluating statistical significance and was
not based on the direction or magnitude of observed group differences.

No additional relative-abundance threshold was applied:

```text
min_abundance = 0
zero_threshold = 0
```

Therefore, a taxon was removed only because of low prevalence and not because
its abundance was below an arbitrary relative-abundance cutoff.

## Normalization and transformation

Taxon counts were normalized independently at each taxonomic rank using
total-sum scaling:

```text
normalization = TSS
```

TSS converts the count profile of each sample into relative abundances and
accounts for differences in the total number of reads when fitting the
abundance component of the model.

Positive relative-abundance values were subsequently log2-transformed:

```text
transform = LOG
```

This is the transformation recommended for the MaAsLin3 abundance model.

The prevalence model was based on whether each taxon was detected or not
detected in each sample.

## Sequencing-depth covariate

The original sequencing depth of each sample was calculated before
normalization as the total number of ASV reads:

```text
library_size = sum of ASV counts
```

It was included in the model as:

```text
log2_library_size = log2(library_size)
```

The variable was standardized by MaAsLin3:

```text
standardize = TRUE
```

Sequencing depth was included because deeper sequencing can increase the
probability of detecting low-prevalence taxa. Without this adjustment,
differences in read depth could be incorrectly interpreted as biological
prevalence differences.

The sequencing-depth covariate was retained in the statistical model even
though it was not displayed in the primary differential-association figures.
Its complete results remain available in the MaAsLin3 result tables.

## Statistical model

The same model was fitted independently at every taxonomic rank:

```r
~ group + log2_library_size
```

The factor order was:

```r
c(
  "hibernation",
  "active",
  "post_hibernation"
)
```

Consequently, MaAsLin3 reported the following group coefficients:

- `group active`: ABH versus H;
- `group post_hibernation`: AAH versus H.

## Abundance and prevalence components

MaAsLin3 evaluates two related components for every taxon.

### Abundance

The abundance model evaluates changes in the amount of a taxon among samples
in which the taxon was observed.

For a categorical group comparison, the abundance coefficient is on a log2
scale. Therefore:

```text
relative-abundance fold change = 2^β
```

### Prevalence

The prevalence model evaluates changes in the probability that a taxon is
detected.

The prevalence coefficient represents a change in log odds. Therefore:

```text
odds ratio = exp(β)
```

Positive coefficients indicate higher abundance or prevalence in ABH or AAH
relative to H. Negative coefficients indicate lower abundance or prevalence
relative to H.

## Compositionality correction

The abundance analysis used:

```text
median_comparison_abundance = TRUE
```

Microbial relative abundances are compositional because the proportions in
each sample sum to a constant. An increase in one taxon can therefore produce
apparent decreases in other taxa even when their absolute abundances have not
changed.

To reduce this problem, abundance coefficients were tested against the median
coefficient across taxa for the same metadata comparison rather than against
zero.

In the custom figures:

- circles represent abundance coefficients;
- the grey dashed vertical line represents the median abundance null
  hypothesis.

The coefficients themselves remain relative-abundance coefficients because
the median was not subtracted from the reported effect sizes.

The prevalence analysis used:

```text
median_comparison_prevalence = FALSE
```

Prevalence coefficients were therefore tested against zero.

In the custom figures:

- triangles represent prevalence coefficients;
- the solid black vertical line represents the prevalence null hypothesis at
  beta equal to zero.

## Prevalence-model stabilization

The following MaAsLin3 options were enabled:

```text
augment = TRUE
warn_prevalence = TRUE
```

Augmentation was used to reduce problems caused by complete or near-complete
separation in logistic prevalence models.

Prevalence warnings were retained so that associations potentially induced by
stronger abundance effects could be identified in the complete result tables.

## Multiple-testing correction

Benjamini-Hochberg false-discovery-rate correction was used:

```text
correction = BH
max_significance = 0.05
```

Each taxonomic rank was analyzed in a separate MaAsLin3 run. Consequently,
multiple-testing correction was performed independently within each
rank-specific analysis rather than across pooled results from all taxonomic
ranks.

The principal significance threshold was:

```text
q < 0.05
```

MaAsLin3 reports two types of corrected results:

- `qval_individual`: significance of the abundance or prevalence component;
- `qval_joint`: joint evidence that either abundance or prevalence is
  associated with the metadata variable.

The complete output retains both individual and joint q-values.

## Visualization strategy

Separate custom figures were produced for phylum, class, order, family, and
genus.

A taxon was included in a figure if at least one group comparison had:

```text
qval_individual < 0.05
```

for either the abundance or prevalence component.

This filtering was applied only to the visualization. It did not remove rows
from the complete MaAsLin3 result tables.

For every included taxon, the figure displays both abundance and prevalence
estimates for both group contrasts. A non-significant component is therefore
still displayed with a white symbol when the other component is significant.

Symbol definitions:

- circle: abundance;
- triangle: prevalence.

Colour definitions:

- purple: significant abundance association;
- green: significant prevalence association;
- white: individual q-value greater than or equal to 0.05.

Confidence intervals represent approximately 95% intervals calculated as:

```text
β ± 1.96 × standard error
```

Taxon names at the family and genus ranks were italicized, while the word
“Unclassified” and unresolved database identifiers were retained in regular
font.

## Output files

The principal combined result files are:

```text
all_ranks_all_results.tsv
all_ranks_group_results.tsv
all_ranks_significant_results.tsv
all_ranks_group_significant_results.tsv
rank_filter_summary.tsv
```

The exact data used for every custom figure are also exported as rank-specific
TSV files.

Generated data, figures, and MaAsLin3 model objects are stored outside the Git
repository. The repository contains only the analysis code and methodological
documentation.

## Reproducibility

The analysis is implemented in:

```text
workflows/16s/r/04_maaslin3_daa.R
workflows/16s/r/05_visualize_maaslin3.R
```

The R and package versions used for model fitting are recorded in:

```text
session_info.txt
```

## References

- MaAsLin3 documentation:
  https://github.com/biobakery/maaslin3
- MaAsLin3 Bioconductor manual:
  https://bioconductor.org/packages/maaslin3
- Nickols WA, Kuntz T, Shen J, et al. MaAsLin 3: refining and extending
  generalized multivariable linear models for meta-omic association discovery.
  Nature Methods. 2026.
  https://doi.org/10.1038/s41592-025-02923-9
