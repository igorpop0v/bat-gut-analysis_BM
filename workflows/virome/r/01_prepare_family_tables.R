# Prepare viral family count tables for downstream analyses.
#
# The main dataset contains all 10 samples without rarefaction.
# RNA_S12046Nr12 is retained but marked as a low-depth sample.
# A second dataset without this sample is created for sensitivity analyses.

library(dplyr)
library(tidyr)
library(readr)

source("workflows/virome/r/local_paths.R")

# Input files.
counts_file <- file.path(
  virome_tables_dir,
  "viral_family_counts_long.tsv"
)

qc_file <- file.path(
  virome_tables_dir,
  "sample_qc.tsv"
)

# Output directory.
output_dir <- file.path(
  virome_results_dir,
  "prepared_tables"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Sample with very low viral sequencing depth.
low_depth_sample <- "RNA_S12046Nr12"

# Experimental group order.
group_levels <- c(
  "active",
  "hibernation",
  "post_hibernation"
)

# Read metadata.
metadata <- read_tsv(
  virome_metadata_file,
  show_col_types = FALSE
) |>
  mutate(
    group = factor(
      group,
      levels = group_levels
    )
  )

# Read sample-level QC information.
sample_qc <- read_tsv(
  qc_file,
  show_col_types = FALSE
)

# Combine metadata and QC.
sample_information <- metadata |>
  left_join(
    sample_qc,
    by = "sample_id"
  ) |>
  mutate(
    low_depth = sample_id == low_depth_sample,
    log2_non_rrna_pairs = log2(metabuli_input_pairs + 1)
  )

write_tsv(
  sample_information |>
    mutate(group = as.character(group)),
  file.path(
    output_dir,
    "sample_metadata_and_qc.tsv"
  )
)

# Read viral family counts.
family_counts_long <- read_tsv(
  counts_file,
  show_col_types = FALSE
) |>
  select(
    sample_id,
    taxid,
    family,
    count
  ) |>
  group_by(
    sample_id,
    family
  ) |>
  summarise(
    count = sum(count),
    .groups = "drop"
  )

# Add explicit zeroes for families not detected in a sample.
family_counts_complete <- expand_grid(
  sample_id = metadata$sample_id,
  family = sort(unique(family_counts_long$family))
) |>
  left_join(
    family_counts_long,
    by = c("sample_id", "family")
  ) |>
  mutate(
    count = replace_na(count, 0)
  )

# Calculate relative abundance within the classified virome.
family_abundance_complete <- family_counts_complete |>
  group_by(sample_id) |>
  mutate(
    viral_family_reads = sum(count),
    relative_abundance = count / viral_family_reads,
    relative_abundance_pct = 100 * relative_abundance
  ) |>
  ungroup()

# Create a family-level summary.
family_summary <- family_abundance_complete |>
  group_by(family) |>
  summarise(
    total_reads = sum(count),
    positive_samples = sum(count > 0),
    prevalence = positive_samples / n_distinct(sample_id),
    mean_relative_abundance_pct = mean(relative_abundance_pct),
    median_relative_abundance_pct = median(relative_abundance_pct),
    maximum_relative_abundance_pct = max(relative_abundance_pct),
    .groups = "drop"
  ) |>
  arrange(
    desc(total_reads)
  )

write_tsv(
  family_summary,
  file.path(
    output_dir,
    "viral_family_summary.tsv"
  )
)

# Complete long-format table.
write_tsv(
  family_abundance_complete,
  file.path(
    output_dir,
    "viral_family_counts_long_complete.tsv"
  )
)

# Count matrix: samples in rows and families in columns.
family_counts_wide <- family_counts_complete |>
  pivot_wider(
    names_from = family,
    values_from = count,
    values_fill = 0
  ) |>
  arrange(
    match(sample_id, metadata$sample_id)
  )

write_tsv(
  family_counts_wide,
  file.path(
    output_dir,
    "viral_family_counts_all_samples.tsv"
  )
)

# Relative-abundance matrix.
family_abundance_wide <- family_abundance_complete |>
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
    match(sample_id, metadata$sample_id)
  )

write_tsv(
  family_abundance_wide,
  file.path(
    output_dir,
    "viral_family_relative_abundance_all_samples.tsv"
  )
)

# Sensitivity dataset without the low-depth sample.
write_tsv(
  family_counts_wide |>
    filter(sample_id != low_depth_sample),
  file.path(
    output_dir,
    "viral_family_counts_without_nr12.tsv"
  )
)

write_tsv(
  family_abundance_wide |>
    filter(sample_id != low_depth_sample),
  file.path(
    output_dir,
    "viral_family_relative_abundance_without_nr12.tsv"
  )
)

write_tsv(
  sample_information |>
    filter(sample_id != low_depth_sample) |>
    mutate(group = as.character(group)),
  file.path(
    output_dir,
    "sample_metadata_without_nr12.tsv"
  )
)

cat("\nViral family tables were prepared.\n")
cat("Rarefaction was not performed.\n")
cat("Main samples:", nrow(sample_information), "\n")
cat(
  "Sensitivity samples:",
  sum(sample_information$sample_id != low_depth_sample),
  "\n"
)
cat("Detected viral families:", nrow(family_summary), "\n")
cat("Results:", output_dir, "\n\n")

print(
  sample_information |>
    select(
      sample_id,
      group_label,
      metabuli_input_pairs,
      viral_pairs,
      low_depth
    )
)
