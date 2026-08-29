# Differential abundance and prevalence analysis of viral families
# using MaAsLin3.
#
# Primary analysis:
#   all samples and all classified viral families.
#
# Sensitivity analyses:
#   1. adjustment for non-rRNA sequencing depth;
#   2. exclusion of the low-depth sample RNA_S12046Nr12;
#   3. exclusion of Parvoviridae;
#   4. exclusion of both RNA_S12046Nr12 and Parvoviridae.
#
# Reference group: hibernation.

library(dplyr)
library(readr)
library(tibble)
library(maaslin3)

source("workflows/virome/r/local_paths.R")

# -------------------------------------------------------------------------
# Input and output paths
# -------------------------------------------------------------------------

prepared_input_dir <- file.path(
  virome_results_dir,
  "prepared_tables"
)

counts_file <- file.path(
  prepared_input_dir,
  "viral_family_counts_all_samples.tsv"
)

metadata_file <- file.path(
  prepared_input_dir,
  "sample_metadata_and_qc.tsv"
)

category_file <- file.path(
  "workflows",
  "virome",
  "config",
  "viral_family_categories.tsv"
)

analysis_root <- file.path(
  virome_results_dir,
  "maaslin3",
  "viral_families"
)

if (dir.exists(analysis_root)) {
  stop(
    "Output directory already exists: ",
    analysis_root,
    "\nRename it before rerunning the analysis."
  )
}

dir.create(
  analysis_root,
  recursive = TRUE,
  showWarnings = FALSE
)

# -------------------------------------------------------------------------
# Analysis settings
# -------------------------------------------------------------------------

low_depth_sample <- "RNA_S12046Nr12"

# A family must have at least 10 reads across the dataset.
min_total_reads <- 10

# A family must occur in at least two samples.
# With only 10 samples, a 10% filter would retain single-sample detections.
min_positive_samples <- 2

# -------------------------------------------------------------------------
# Read input data
# -------------------------------------------------------------------------

family_counts <- read_tsv(
  counts_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

metadata <- read_tsv(
  metadata_file,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
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

family_categories <- read_tsv(
  category_file,
  show_col_types = FALSE
)

if (any(is.na(metadata$group))) {
  stop("Unknown or missing group labels were found.")
}

if (!setequal(
  family_counts$sample_id,
  metadata$sample_id
)) {
  stop("Sample IDs differ between the count table and metadata.")
}

if (!"log2_non_rrna_pairs" %in% names(metadata)) {
  stop(
    "The metadata table does not contain log2_non_rrna_pairs."
  )
}

# -------------------------------------------------------------------------
# Models to run
# -------------------------------------------------------------------------

model_plan <- tribble(
  ~analysis, ~remove_nr12, ~remove_parvoviridae, ~adjust_for_depth,
  "primary_all_samples", FALSE, FALSE, FALSE,
  "depth_adjusted", FALSE, FALSE, TRUE,
  "without_nr12", TRUE, FALSE, FALSE,
  "without_parvoviridae", FALSE, TRUE, FALSE,
  "without_nr12_and_parvoviridae", TRUE, TRUE, FALSE
)

# -------------------------------------------------------------------------
# Function for running one MaAsLin3 model
# -------------------------------------------------------------------------

run_virome_model <- function(
    analysis_name,
    remove_nr12,
    remove_parvoviridae,
    adjust_for_depth) {

  message("")
  message("Running analysis: ", analysis_name)

  analysis_dir <- file.path(
    analysis_root,
    analysis_name
  )

  prepared_dir <- file.path(
    analysis_dir,
    "prepared_inputs"
  )

  model_dir <- file.path(
    analysis_dir,
    "model"
  )

  dir.create(
    prepared_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Select samples for this model.
  current_metadata <- metadata

  if (remove_nr12) {
    current_metadata <- current_metadata |>
      filter(sample_id != low_depth_sample)
  }

  current_metadata <- current_metadata |>
    arrange(
      match(
        sample_id,
        metadata$sample_id
      )
    )

  current_counts <- family_counts |>
    filter(
      sample_id %in% current_metadata$sample_id
    ) |>
    arrange(
      match(
        sample_id,
        current_metadata$sample_id
      )
    )

  # Remove the non-taxonomic category.
  current_counts <- current_counts |>
    select(
      -any_of("Unclassified viruses")
    )

  # Remove Parvoviridae only in the specified
  # sensitivity analyses.
  if (remove_parvoviridae) {
    current_counts <- current_counts |>
      select(
        -any_of("Parvoviridae")
      )
  }

  count_matrix <- as.matrix(
    current_counts[, -1, drop = FALSE]
  )

  storage.mode(count_matrix) <- "numeric"

  if (any(count_matrix < 0)) {
    stop(
      "Negative counts found in analysis: ",
      analysis_name
    )
  }

  # Calculate prevalence and total read support.
  filter_audit <- tibble(
    family = colnames(count_matrix),
    total_reads = colSums(count_matrix),
    positive_samples = colSums(count_matrix > 0),
    total_samples = nrow(count_matrix),
    prevalence = positive_samples / total_samples
  ) |>
    mutate(
      passes_filter =
        total_reads >= min_total_reads &
        positive_samples >= min_positive_samples
    ) |>
    left_join(
      family_categories,
      by = "family"
    ) |>
    arrange(
      desc(passes_filter),
      desc(total_reads)
    )

  retained_families <- filter_audit |>
    filter(passes_filter) |>
    pull(family)

  if (length(retained_families) == 0) {
    stop(
      "No viral families passed filtering in analysis: ",
      analysis_name
    )
  }

  filtered_matrix <- count_matrix[
    ,
    retained_families,
    drop = FALSE
  ]

  feature_table <- as.data.frame(
    filtered_matrix,
    check.names = FALSE
  )

  rownames(feature_table) <- current_counts$sample_id

  # Primary formula contains only the biological group.
  # Sequencing depth is included only in a sensitivity model.
  if (adjust_for_depth) {
    model_formula <- ~ group + log2_non_rrna_pairs
    model_formula_text <- "~ group + log2_non_rrna_pairs"
  } else {
    model_formula <- ~ group
    model_formula_text <- "~ group"
  }

  metadata_columns <- c(
    "sample_id",
    "group"
  )

  if (adjust_for_depth) {
    metadata_columns <- c(
      metadata_columns,
      "log2_non_rrna_pairs"
    )
  }

  metadata_table <- current_metadata |>
    select(
      all_of(metadata_columns)
    ) |>
    column_to_rownames("sample_id") |>
    as.data.frame()

  if (!identical(
    rownames(feature_table),
    rownames(metadata_table)
  )) {
    stop(
      "Sample order differs in analysis: ",
      analysis_name
    )
  }

  # Save exactly what enters MaAsLin3.
  write_tsv(
    feature_table |>
      rownames_to_column("sample_id"),
    file.path(
      prepared_dir,
      "viral_family_counts_filtered.tsv"
    )
  )

  write_tsv(
    current_metadata |>
      mutate(
        group = as.character(group)
      ) |>
      select(
        all_of(metadata_columns)
      ),
    file.path(
      prepared_dir,
      "sample_metadata_maaslin3.tsv"
    )
  )

  write_tsv(
    filter_audit,
    file.path(
      prepared_dir,
      "viral_family_filter_audit.tsv"
    )
  )

  settings <- tibble(
    parameter = c(
      "Analysis",
      "Samples",
      "Families before filtering",
      "Families after filtering",
      "Minimum total reads",
      "Minimum positive samples",
      "Normalization",
      "Transformation",
      "Reference group",
      "Model formula",
      "RNA_S12046Nr12 excluded",
      "Parvoviridae excluded"
    ),
    value = c(
      analysis_name,
      nrow(feature_table),
      ncol(count_matrix),
      ncol(feature_table),
      min_total_reads,
      min_positive_samples,
      "TSS",
      "LOG",
      "hibernation",
      model_formula_text,
      remove_nr12,
      remove_parvoviridae
    )
  )

  write_tsv(
    settings,
    file.path(
      analysis_dir,
      "analysis_settings.tsv"
    )
  )

  # -----------------------------------------------------------------------
  # Run MaAsLin3
  # -----------------------------------------------------------------------

  set.seed(12345)

  maaslin3(
    input_data = feature_table,
    input_metadata = metadata_table,
    output = model_dir,

    formula = model_formula,

    # Relative abundance normalization.
    normalization = "TSS",

    # Log2 transformation of positive relative abundances.
    transform = "LOG",

    # Filtering was performed explicitly above.
    min_prevalence = 0,
    min_abundance = 0,
    zero_threshold = 0,

    # Account for compositionality in the abundance model.
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE,

    # Evaluate both abundance and prevalence.
    warn_prevalence = TRUE,
    augment = TRUE,

    # Benjamini-Hochberg correction.
    correction = "BH",
    max_significance = 0.05,

    standardize = TRUE,

    # Standard MaAsLin3 output for initial inspection.
    plot_summary_plot = TRUE,
    summary_plot_first_n = 25,
    plot_associations = FALSE,

    cores = 1
  )

  # -----------------------------------------------------------------------
  # Prepare readable result tables
  # -----------------------------------------------------------------------

  all_results <- read_tsv(
    file.path(
      model_dir,
      "all_results.tsv"
    ),
    show_col_types = FALSE
  ) |>
    left_join(
      family_categories,
      by = c("feature" = "family")
    ) |>
    mutate(
      analysis = analysis_name,
      .before = 1
    )

  group_results <- all_results |>
    filter(metadata == "group") |>
    mutate(
      comparison = case_when(
        value == "active" ~ "ABH vs H",
        value == "post_hibernation" ~ "AAH vs H",
        TRUE ~ as.character(value)
      ),
      higher_in = case_when(
        comparison == "ABH vs H" & coef > 0 ~ "ABH",
        comparison == "ABH vs H" & coef < 0 ~ "H",
        comparison == "AAH vs H" & coef > 0 ~ "AAH",
        comparison == "AAH vs H" & coef < 0 ~ "H",
        TRUE ~ NA_character_
      )
    ) |>
    relocate(
      comparison,
      higher_in,
      .after = value
    )

  significant_group_results <- group_results |>
    filter(
      !is.na(qval_individual),
      qval_individual < 0.05
    ) |>
    arrange(
      comparison,
      model,
      qval_individual
    )

  write_tsv(
    all_results,
    file.path(
      analysis_dir,
      "all_results_annotated.tsv"
    )
  )

  write_tsv(
    group_results,
    file.path(
      analysis_dir,
      "group_results.tsv"
    )
  )

  write_tsv(
    significant_group_results,
    file.path(
      analysis_dir,
      "group_significant_results_q05.tsv"
    )
  )

  model_summary <- tibble(
    analysis = analysis_name,
    samples = nrow(feature_table),
    ABH_samples = sum(
      current_metadata$group == "active"
    ),
    H_samples = sum(
      current_metadata$group == "hibernation"
    ),
    AAH_samples = sum(
      current_metadata$group == "post_hibernation"
    ),
    families_before_filtering = ncol(count_matrix),
    families_after_filtering = ncol(feature_table),
    significant_group_associations =
      nrow(significant_group_results),
    significant_families =
      n_distinct(significant_group_results$feature),
    formula = model_formula_text
  )

  return(
    list(
      summary = model_summary,
      all_results = all_results,
      group_results = group_results,
      significant_results = significant_group_results
    )
  )
}

# -------------------------------------------------------------------------
# Run all models
# -------------------------------------------------------------------------

model_results <- vector(
  "list",
  nrow(model_plan)
)

for (i in seq_len(nrow(model_plan))) {

  model_results[[i]] <- run_virome_model(
    analysis_name = model_plan$analysis[[i]],
    remove_nr12 = model_plan$remove_nr12[[i]],
    remove_parvoviridae =
      model_plan$remove_parvoviridae[[i]],
    adjust_for_depth =
      model_plan$adjust_for_depth[[i]]
  )
}

# -------------------------------------------------------------------------
# Combine results from all models
# -------------------------------------------------------------------------

model_summary <- bind_rows(
  lapply(
    model_results,
    function(x) x$summary
  )
)

combined_group_results <- bind_rows(
  lapply(
    model_results,
    function(x) x$group_results
  )
)

combined_significant_results <- bind_rows(
  lapply(
    model_results,
    function(x) x$significant_results
  )
)

write_tsv(
  model_summary,
  file.path(
    analysis_root,
    "model_summary.tsv"
  )
)

write_tsv(
  combined_group_results,
  file.path(
    analysis_root,
    "all_models_group_results.tsv"
  )
)

write_tsv(
  combined_significant_results,
  file.path(
    analysis_root,
    "all_models_group_significant_results_q05.tsv"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    analysis_root,
    "session_info.txt"
  )
)

primary_significant <- combined_significant_results |>
  filter(
    analysis == "primary_all_samples"
  )

cat(
  "\nVirome MaAsLin3 analysis completed.\n",
  "Primary significant group associations: ",
  nrow(primary_significant),
  "\nResults: ",
  analysis_root,
  "\n",
  sep = ""
)

print(model_summary)
