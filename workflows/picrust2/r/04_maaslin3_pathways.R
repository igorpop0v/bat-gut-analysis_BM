# Differential abundance and prevalence analysis of predicted MetaCyc pathways.
# Reference group: hibernation.
#
# Comparisons:
#   active vs hibernation
#   post_hibernation vs hibernation
#
# MaAsLin3 evaluates:
#   1. abundance among samples where a pathway is detected;
#   2. prevalence, meaning the probability that a pathway is detected.

library(dplyr)
library(readr)
library(tibble)
library(maaslin3)

source("workflows/picrust2/r/local_paths.R")

# Input files.
pathway_file <- file.path(
  picrust2_input_dir,
  "metacyc_pathway_abundance_with_descriptions.tsv.gz"
)

metadata_file <- file.path(
  picrust2_input_dir,
  "sample_metadata.tsv"
)

asv_counts_file <- file.path(
  picrust2_input_dir,
  "asv_counts_unrarefied.tsv"
)

# Output directories.
analysis_root <- file.path(
  picrust2_results_dir,
  "maaslin3",
  "pathways_primary"
)

prepared_dir <- file.path(
  analysis_root,
  "prepared_inputs"
)

model_output_dir <- file.path(
  analysis_root,
  "model"
)

# Prevent accidental overwriting of existing results.
if (dir.exists(analysis_root)) {
  stop(
    "Output directory already exists: ",
    analysis_root,
    "\nRename it before rerunning the analysis."
  )
}

dir.create(
  prepared_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# A pathway must be detected in at least 10% of samples.
min_prevalence <- 0.10

# -------------------------------------------------------------------------
# 1. Read predicted MetaCyc pathway abundances
# -------------------------------------------------------------------------

pathway_input <- read_tsv(
  pathway_file,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  mutate(
    description = if_else(
      is.na(description) | description == "not_found",
      pathway,
      description
    )
  )

if (anyDuplicated(pathway_input$pathway)) {
  stop("Duplicated MetaCyc pathway IDs were found.")
}

pathway_descriptions <- pathway_input |>
  select(pathway, description)

write_tsv(
  pathway_descriptions,
  file.path(
    prepared_dir,
    "metacyc_pathway_descriptions.tsv"
  )
)

# -------------------------------------------------------------------------
# 2. Read sample metadata
# -------------------------------------------------------------------------

metadata <- read_tsv(
  metadata_file,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  rename(sample_id = 1) |>
  mutate(
    group = factor(
      group,
      levels = c(
        "hibernation",
        "active",
        "post_hibernation"
      )
    )
  )

if (any(is.na(metadata$group))) {
  stop("Unknown or missing group labels were found in the metadata.")
}

if (!all(metadata$sample_id %in% names(pathway_input))) {
  missing_samples <- setdiff(
    metadata$sample_id,
    names(pathway_input)
  )

  stop(
    "Samples missing from the pathway table: ",
    paste(missing_samples, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# 3. Create a samples-by-pathways feature table
# -------------------------------------------------------------------------

feature_table <- t(
  as.matrix(
    pathway_input[, metadata$sample_id]
  )
) |>
  as.data.frame(check.names = FALSE)

rownames(feature_table) <- metadata$sample_id
colnames(feature_table) <- pathway_input$pathway

if (any(as.matrix(feature_table) < 0)) {
  stop("Negative predicted pathway abundances were found.")
}

write_tsv(
  feature_table |>
    rownames_to_column("sample_id"),
  file.path(
    prepared_dir,
    "pathway_abundance_unnormalized.tsv"
  )
)

# -------------------------------------------------------------------------
# 4. Filter pathways by prevalence
# -------------------------------------------------------------------------

feature_matrix <- as.matrix(feature_table)

positive_samples <- colSums(feature_matrix > 0)

prevalence_table <- tibble(
  pathway = colnames(feature_matrix),
  positive_samples = positive_samples,
  total_samples = nrow(feature_matrix),
  prevalence = positive_samples / nrow(feature_matrix)
) |>
  left_join(
    pathway_descriptions,
    by = "pathway"
  ) |>
  mutate(
    passes_filter = prevalence >= min_prevalence
  ) |>
  arrange(
    desc(prevalence),
    pathway
  )

write_tsv(
  prevalence_table,
  file.path(
    prepared_dir,
    "pathway_prevalence_filter_audit.tsv"
  )
)

retained_pathways <- prevalence_table |>
  filter(passes_filter) |>
  pull(pathway)

feature_table_filtered <- feature_table[
  ,
  retained_pathways,
  drop = FALSE
]

write_tsv(
  feature_table_filtered |>
    rownames_to_column("sample_id"),
  file.path(
    prepared_dir,
    "pathway_abundance_after_prevalence_filter.tsv"
  )
)

# -------------------------------------------------------------------------
# 5. Calculate original 16S sequencing depth
# -------------------------------------------------------------------------

asv_counts <- read_tsv(
  asv_counts_file,
  skip = 1,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  rename(feature_id = 1)

if (!all(metadata$sample_id %in% names(asv_counts))) {
  missing_samples <- setdiff(
    metadata$sample_id,
    names(asv_counts)
  )

  stop(
    "Samples missing from the ASV count table: ",
    paste(missing_samples, collapse = ", ")
  )
}

library_sizes <- colSums(
  as.matrix(
    asv_counts[, metadata$sample_id]
  )
)

model_metadata <- metadata |>
  mutate(
    library_size = as.numeric(
      library_sizes[sample_id]
    ),
    log2_library_size = log2(library_size)
  )

if (any(!is.finite(model_metadata$log2_library_size))) {
  stop("Invalid 16S library sizes were detected.")
}

metadata_table <- model_metadata |>
  select(
    sample_id,
    group,
    log2_library_size
  ) |>
  column_to_rownames("sample_id") |>
  as.data.frame()

write_tsv(
  model_metadata |>
    select(
      sample_id,
      group,
      library_size,
      log2_library_size
    ),
  file.path(
    prepared_dir,
    "maaslin3_metadata.tsv"
  )
)

# -------------------------------------------------------------------------
# 6. Save a filtering and modelling summary
# -------------------------------------------------------------------------

minimum_positive_samples <- ceiling(
  min_prevalence * nrow(feature_table)
)

filter_summary <- tibble(
  parameter = c(
    "Samples",
    "Pathways before filtering",
    "Pathways after filtering",
    "Pathways removed",
    "Minimum prevalence",
    "Minimum positive samples",
    "Normalization",
    "Transformation",
    "Reference group",
    "Model formula"
  ),
  value = c(
    nrow(feature_table),
    ncol(feature_table),
    ncol(feature_table_filtered),
    ncol(feature_table) -
      ncol(feature_table_filtered),
    min_prevalence,
    minimum_positive_samples,
    "TSS",
    "LOG",
    "hibernation",
    "~ group + log2_library_size"
  )
)

write_tsv(
  filter_summary,
  file.path(
    analysis_root,
    "pathway_filter_summary.tsv"
  )
)

# -------------------------------------------------------------------------
# 7. Run MaAsLin3
# -------------------------------------------------------------------------

set.seed(12345)

fit <- maaslin3(
  input_data = feature_table_filtered,
  input_metadata = metadata_table,
  output = model_output_dir,

  formula = ~ group + log2_library_size,

  normalization = "TSS",
  transform = "LOG",

  min_prevalence = 0,
  min_abundance = 0,
  zero_threshold = 0,

  median_comparison_abundance = TRUE,
  median_comparison_prevalence = FALSE,

  warn_prevalence = TRUE,
  augment = TRUE,

  correction = "BH",
  max_significance = 0.05,
  standardize = TRUE,

  plot_summary_plot = TRUE,
  plot_associations = FALSE,

  cores = 1
)

# -------------------------------------------------------------------------
# 8. Add pathway descriptions to MaAsLin3 results
# -------------------------------------------------------------------------

annotate_results <- function(file_name) {

  read_tsv(
    file.path(model_output_dir, file_name),
    show_col_types = FALSE
  ) |>
    left_join(
      pathway_descriptions,
      by = c("feature" = "pathway")
    ) |>
    relocate(
      description,
      .after = feature
    )
}

all_results <- annotate_results(
  "all_results.tsv"
)

software_significant_results <- annotate_results(
  "significant_results.tsv"
)

# Keep group comparisons and give them clear labels.
group_results <- all_results |>
  filter(metadata == "group") |>
  mutate(
    comparison = case_when(
      value == "active" ~ "ABH vs H",
      value == "post_hibernation" ~ "AAH vs H",
      TRUE ~ as.character(value)
    )
  ) |>
  relocate(
    comparison,
    .after = value
  )

# Primary significant results:
# individual abundance or prevalence q-value below 0.05.
group_significant_results <- group_results |>
  filter(qval_individual < 0.05) |>
  arrange(
    comparison,
    model,
    qval_individual
  )

write_tsv(
  all_results,
  file.path(
    analysis_root,
    "pathway_all_results_with_descriptions.tsv"
  )
)

write_tsv(
  software_significant_results,
  file.path(
    analysis_root,
    "pathway_software_significant_results_with_descriptions.tsv"
  )
)

write_tsv(
  group_results,
  file.path(
    analysis_root,
    "pathway_group_results_with_descriptions.tsv"
  )
)

write_tsv(
  group_significant_results,
  file.path(
    analysis_root,
    "pathway_group_significant_results_q05.tsv"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    analysis_root,
    "session_info.txt"
  )
)

cat(
  "\nPICRUSt2 pathway MaAsLin3 analysis completed.\n",
  "Samples: ", nrow(feature_table_filtered), "\n",
  "Pathways before filtering: ", ncol(feature_table), "\n",
  "Pathways after filtering: ", ncol(feature_table_filtered), "\n",
  "Significant group associations: ",
  nrow(group_significant_results), "\n",
  "Results: ", analysis_root, "\n",
  sep = ""
)
