# Virome diversity and sensitivity analysis

## Purpose

This analysis evaluates alpha and beta diversity of the viral-family composition detected by MetaBuli in faecal RNA-seq libraries from common noctules (*Nyctalus noctula*).

The study includes three independent groups:

- **ABH** — active bats before hibernation (`n = 4`);
- **H** — bats during hibernation (`n = 3`);
- **AAH** — active bats after hibernation (`n = 3`).

Because the virome dataset contains only ten samples and sequencing depth differs substantially among libraries, all statistical results are treated as exploratory.

## Input data

The analysis uses viral-family read counts obtained after:

1. read quality control and adapter trimming;
2. removal of reads mapping to the host reference;
3. removal of ribosomal RNA reads with SortMeRNA;
4. taxonomic classification of the remaining paired reads with MetaBuli.

`Unclassified viruses` were excluded from the diversity analysis because this category does not represent a defined viral family.

Families supported by fewer than 10 reads across the entire dataset were excluded to reduce the influence of isolated and potentially unreliable classifications.

## Sequencing-depth limitation

The number of virus-classified read pairs varied substantially among samples. In particular, sample `RNA_S12046Nr12` contained only 164 viral read pairs.

Rarefying all samples to 164 reads would discard nearly all information from the more deeply sequenced libraries. It would also introduce additional random subsampling variation while retaining very limited taxonomic information.

Therefore, rarefaction was not used for the primary virome diversity analysis. Instead:

- alpha-diversity indices were calculated from within-sample relative abundances;
- relative abundances were Hellinger transformed before calculation of Bray–Curtis dissimilarities;
- the robustness of the results was evaluated through sensitivity analyses.

This decision applies only to the exploratory virome dataset and should not be interpreted as a general recommendation against rarefaction for every microbiome analysis.

## Primary analysis

The primary analysis retained:

- all ten samples, including the low-depth sample `RNA_S12046Nr12`;
- all supported viral families, including *Parvoviridae*.

The low-depth sample was retained to avoid post hoc sample exclusion. *Parvoviridae* was retained because it represents a genuine and biologically prominent component of the detected viral composition.

The following alpha-diversity indices were calculated:

- Shannon index;
- Simpson index;
- Pielou evenness.

Differences among the three groups were tested using the Kruskal–Wallis test. Pairwise Wilcoxon tests were adjusted using the Benjamini–Hochberg procedure.

Beta diversity was evaluated using:

- total-sum scaling to relative abundance;
- Hellinger transformation;
- Bray–Curtis dissimilarity;
- principal coordinates analysis;
- PERMANOVA with 9,999 permutations;
- PERMDISP with 9,999 permutations.

## Sensitivity analyses

Four versions of the dataset were compared:

1. **Primary dataset** — all samples and all supported viral families;
2. **Without the low-depth sample** — `RNA_S12046Nr12` excluded;
3. **Without Parvoviridae** — all samples retained, but *Parvoviridae* excluded;
4. **Without both** — `RNA_S12046Nr12` and *Parvoviridae* excluded.

The exclusion of *Parvoviridae* is a sensitivity analysis rather than part of the primary filtering strategy. This family may contain a substantial insect-associated signal in the faecal samples of insectivorous bats and strongly dominates several libraries.

At the family-classification level, the complete *Parvoviridae* signal should not automatically be described as densovirus. Evidence for densovirus identity must instead be supported by sequence assembly and genome-level analysis.

## Sensitivity-analysis results

| Dataset | Samples | Shannon p | Simpson p | Pielou p | PERMANOVA R² | PERMANOVA p | PERMDISP p |
|---|---:|---:|---:|---:|---:|---:|---:|
| Primary dataset | 10 | 0.043 | 0.043 | 0.055 | 0.292 | 0.156 | 0.503 |
| Without `RNA_S12046Nr12` | 9 | 0.086 | 0.086 | 0.118 | 0.341 | 0.108 | 0.438 |
| Without *Parvoviridae* | 10 | 0.921 | 0.705 | 0.689 | 0.302 | 0.118 | 0.981 |
| Without both | 9 | 0.962 | 0.951 | 0.741 | 0.293 | 0.225 | 0.965 |

In the primary analysis, Shannon and Simpson indices produced nominal global p-values below 0.05. However:

- none of the pairwise comparisons remained significant after multiple-testing correction;
- the global alpha-diversity results were no longer significant after excluding the low-depth sample;
- the apparent differences disappeared after excluding *Parvoviridae*.

Consequently, the alpha-diversity result is not considered robust. It appears to be strongly influenced by the dominance of *Parvoviridae* and, to a lesser extent, by the low-depth sample.

The global PERMANOVA result was not significant in the primary analysis or in any sensitivity analysis. All adjusted pairwise PERMANOVA comparisons were also non-significant.

PERMDISP was non-significant in every analysis, providing no evidence that differences in within-group dispersion explain the ordination pattern.

## Interpretation

The analysis does not provide robust evidence that hibernation status is associated with overall viral-family alpha or beta diversity.

The nominal alpha-diversity result in the complete dataset should be interpreted cautiously because it is not preserved across sensitivity analyses. It should not be reported as a definitive difference among physiological groups.

The most defensible interpretation is that:

- the viral composition is highly heterogeneous among individual samples;
- several samples are dominated by *Parvoviridae*;
- sequencing depth and viral read recovery vary substantially;
- the sample size is too small for confirmatory community-level inference.

The results remain useful as an exploratory description of the faecal RNA virome and for identifying viral targets for subsequent genome-level investigation.

## Phage subset

Phage-family counts were also generated separately. However, phage read recovery was extremely uneven:

- some samples contained only a few phage-classified reads;
- `RNA_S12046Nr12` contained no phage-classified reads;
- the hibernation group was largely represented by a single phage-rich sample.

Therefore, phage composition should be presented descriptively. Formal diversity testing or differential-abundance modelling of the phage-only subset would not be statistically reliable with the current data.

## Reproducibility

The primary diversity analysis is implemented in:

`workflows/virome/r/04_virome_diversity_primary.R`

The sensitivity analysis is implemented in:

`workflows/virome/r/05_virome_diversity_sensitivity.R`

## Reference

McMurdie PJ, Holmes S. Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Computational Biology*. 2014;10(4):e1003531. https://doi.org/10.1371/journal.pcbi.1003531
