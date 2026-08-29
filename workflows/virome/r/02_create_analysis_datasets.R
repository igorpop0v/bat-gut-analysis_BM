# Create separate virome datasets for downstream analyses:
#
# 1. Bacteriophage-only dataset.
# 2. Whole virome excluding Parvoviridae.
# 3. Corresponding sensitivity datasets without RNA_S12046Nr12.
#
# Relative abundances are recalculated separately within each dataset.
# Rarefaction is not performed.

library(dplyr)
library(tidyr)
library(readr)

source("workflows/virome/r/local_paths.R")

prepared_dir <- file.path(
  virome_results_dir,
  "prepared_tables"
)

output_dir <- file.path(
  virome_results_dir,
  "analysis_datasets"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

low_depth_sample <- "RNA_S12046Nr12"

counts_file <- file.path(
  prepared_dir,
  "viral_family_counts_long_complete.tsv"
)

sample_information_file <- file.path(
  prepared_dir,
  "sample_metadata_and_qc.tsv"
)

category_file <- paste0(
  "workflows/virome/config/",
  "viral_family_categories.tsv"
)

counts_long <- read_tsv(
  counts_file,
  show_col_types = FALSE
)

sample_information <- read_tsv(
  sample_information_file,
  show_col_types = FALSE
)

family_categories <- read_tsv(
  category_file,
  show_col_types = FALSE
)

# Add the manually curated host category.
classified_counts <- counts_long |>
  select(
    sample_id,
    family,
    count
  ) |>
  left_join(
    family_categories,
    by = "family"
  )

# Stop if a newly detected family has not been categorized.
missing_categories <- classified_counts |>
  filter(is.na(category)) |>
  distinct(family)

if (nrow(missing_categories) > 0) {
  print(missing_categories)

  stop(
    "Some viral families are missing from viral_family_categories.tsv."
  )
}

write_tsv(
  classified_counts,
  file.path(
    output_dir,
    "viral_family_counts_with_categories.tsv"
  )
)

# ============================================================
# Bacteriophage-only dataset
# ============================================================

phage_counts_long <- classified_counts |>
  filter(category == "bacteriophage") |>
  group_by(sample_id) |>
  mutate(
    phage_reads = sum(count),
    relative_abundance = if_else(
      phage_reads > 0,
      count / phage_reads,
      0
    )
  ) |>
  ungroup()

phage_counts_wide <- phage_counts_long |>
  select(
    sample_id,
    family,
    count
  ) |>
  pivot_wider(
    names_from = family,
    values_from = count,
    values_fill = 0
  ) |>
  arrange(
    match(
      sample_id,
      sample_information$sample_id
    )
  )

phage_abundance_wide <- phage_counts_long |>
  select(
    sample_id,
    family,
    relative_abundance
  ) |>
  pivot_wider(
    names_from = family,
    values_from = relative_abundance,
    values_fill = 0
  ) |>
  arrange(
    match(
      sample_id,
      sample_information$sample_id
    )
  )

write_tsv(
  phage_counts_wide,
  file.path(
    output_dir,
    "phage_family_counts_all_samples.tsv"
  )
)

write_tsv(
  phage_abundance_wide,
  file.path(
    output_dir,
    "phage_family_relative_abundance_all_samples.tsv"
  )
)

write_tsv(
  phage_counts_wide |>
    filter(sample_id != low_depth_sample),
  file.path(
    output_dir,
    "phage_family_counts_without_nr12.tsv"
  )
)

write_tsv(
  phage_abundance_wide |>
    filter(sample_id != low_depth_sample),
  file.path(
    output_dir,
    "phage_family_relative_abundance_without_nr12.tsv"
  )
)

# Summarize bacteriophage sequencing depth.
phage_sample_summary <- phage_counts_long |>
  group_by(sample_id) |>
  summarise(
    phage_reads = sum(count),
    detected_phage_families = sum(count > 0),
    .groups = "drop"
  ) |>
  right_join(
    sample_information,
    by = "sample_id"
  ) |>
  mutate(
    phage_reads = replace_na(phage_reads, 0),
    detected_phage_families = replace_na(
      detected_phage_families,
      0
    ),
    phage_pct_of_viral_reads = if_else(
      viral_pairs > 0,
      100 * phage_reads / viral_pairs,
      0
    )
  ) |>
  select(
    sample_id,
    group,
    group_label,
    viral_pairs,
    phage_reads,
    phage_pct_of_viral_reads,
    detected_phage_families,
    low_depth
  )

write_tsv(
  phage_sample_summary,
  file.path(
    output_dir,
    "phage_sample_summary.tsv"
  )
)

# ============================================================
# Whole virome excluding Parvoviridae
# ============================================================

non_parvovirus_long <- classified_counts |>
  filter(family != "Parvoviridae") |>
  group_by(sample_id) |>
  mutate(
    non_parvovirus_reads = sum(count),
    relative_abundance = if_else(
      non_parvovirus_reads > 0,
      count / non_parvovirus_reads,
      0
    )
  ) |>
  ungroup()

non_parvovirus_counts_wide <- non_parvovirus_long |>
  select(
    sample_id,
    family,
    count
  ) |>
  pivot_wider(
    names_from = family,
    values_from = count,
    values_fill = 0
  ) |>
  arrange(
    match(
      sample_id,
      sample_information$sample_id
    )
  )

non_parvovirus_abundance_wide <- non_parvovirus_long |>
  select(
    sample_id,
    family,
    relative_abundance
  ) |>
  pivot_wider(
    names_from = family,
    values_from = relative_abundance,
    values_fill = 0
  ) |>
  arrange(
    match(
      sample_id,
      sample_information$sample_id
    )
  )

write_tsv(
  non_parvovirus_counts_wide,
  file.path(
    output_dir,
    "viral_family_counts_without_parvoviridae.tsv"
  )
)

write_tsv(
  non_parvovirus_abundance_wide,
  file.path(
    output_dir,
    "viral_family_relative_abundance_without_parvoviridae.tsv"
  )
)

write_tsv(
  non_parvovirus_counts_wide |>
    filter(sample_id != low_depth_sample),
  file.path(
    output_dir,
    "viral_family_counts_without_parvoviridae_and_nr12.tsv"
  )
)

write_tsv(
  non_parvovirus_abundance_wide |>
    filter(sample_id != low_depth_sample),
  file.path(
    output_dir,
    "viral_family_relative_abundance_without_parvoviridae_and_nr12.tsv"
  )
)

cat("\nAnalysis datasets were created.\n")
cat("Rarefaction was not performed.\n")
cat("Results:", output_dir, "\n\n")

print(phage_sample_summary)
