# Method-sensitivity analysis of viral-family associations
# using MaAsLin2.
#
# Reference group: hibernation.

library(dplyr)
library(readr)
library(tibble)
library(Maaslin2)

source("workflows/virome/r/local_paths.R")

# -------------------------------------------------------------------------
# Input and output directories
# -------------------------------------------------------------------------

maaslin3_root <- file.path(
  virome_results_dir,
  "maaslin3",
  "viral_families"
)

analysis_root <- file.path(
  virome_results_dir,
  "maaslin2",
  "viral_families_method_sensitivity"
)

if (!dir.exists(maaslin3_root)) {
  stop(
    "MaAsLin3 results directory was not found: ",
    maaslin3_root
  )
}

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

category_file <- file.path(
  "workflows",
  "virome",
  "config",
  "viral_family_categories.tsv"
)

family_categories <- read_tsv(
  category_file,
  show_col_types = FALSE
)

# -------------------------------------------------------------------------
# Use the same five datasets as in MaAsLin3
# -------------------------------------------------------------------------

model_plan <- tibble(
  analysis = c(
    "primary_all_samples",
    "depth_adjusted",
    "without_nr12",
    "without_parvoviridae",
    "without_nr12_and_parvoviridae"
  ),
  adjust_for_depth = c(
    FALSE,
    TRUE,
    FALSE,
    FALSE,
    FALSE
  )
)

model_summaries <- list()
all_model_results <- list()
all_group_results <- list()
all_significant_results <- list()

# -------------------------------------------------------------------------
# Run MaAsLin2
# -------------------------------------------------------------------------

for (i in seq_len(nrow(model_plan))) {

  analysis_name <- model_plan$analysis[[i]]
  adjust_for_depth <- model_plan$adjust_for_depth[[i]]

  message("")
  message("Running MaAsLin2 analysis: ", analysis_name)

  maaslin3_prepared_dir <- file.path(
    maaslin3_root,
    analysis_name,
    "prepared_inputs"
  )

  counts_file <- file.path(
    maaslin3_prepared_dir,
    "viral_family_counts_filtered.tsv"
  )

  metadata_file <- file.path(
    maaslin3_prepared_dir,
    "sample_metadata_maaslin3.tsv"
  )

  if (!file.exists(counts_file)) {
    stop("Count table was not found: ", counts_file)
  }

  if (!file.exists(metadata_file)) {
    stop("Metadata file was not found: ", metadata_file)
  }

  # Read the exact filtered table used by MaAsLin3.
  counts <- read_tsv(
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

  if (any(is.na(metadata$group))) {
    stop(
      "Unknown group labels in analysis: ",
      analysis_name
    )
  }

  if (!setequal(
    counts$sample_id,
    metadata$sample_id
  )) {
    stop(
      "Sample IDs differ in analysis: ",
      analysis_name
    )
  }

  metadata <- metadata |>
    arrange(
      match(
        sample_id,
        counts$sample_id
      )
    )

  feature_table <- counts |>
    column_to_rownames("sample_id") |>
    as.data.frame(check.names = FALSE)

  if (adjust_for_depth) {

    if (!"log2_non_rrna_pairs" %in% names(metadata)) {
      stop(
        "log2_non_rrna_pairs is missing from analysis: ",
        analysis_name
      )
    }

    metadata_table <- metadata |>
      select(
        sample_id,
        group,
        log2_non_rrna_pairs
      ) |>
      column_to_rownames("sample_id") |>
      as.data.frame()

    fixed_effects <- c(
      "group",
      "log2_non_rrna_pairs"
    )

    model_formula <- "~ group + log2_non_rrna_pairs"

  } else {

    metadata_table <- metadata |>
      select(
        sample_id,
        group
      ) |>
      column_to_rownames("sample_id") |>
      as.data.frame()

    fixed_effects <- "group"

    model_formula <- "~ group"
  }

  if (!identical(
    rownames(feature_table),
    rownames(metadata_table)
  )) {
    stop(
      "Sample order differs in analysis: ",
      analysis_name
    )
  }

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

  # Save an audit copy of the exact inputs.
  write_tsv(
    feature_table |>
      rownames_to_column("sample_id"),
    file.path(
      prepared_dir,
      "viral_family_counts_filtered.tsv"
    )
  )

  write_tsv(
    metadata_table |>
      rownames_to_column("sample_id"),
    file.path(
      prepared_dir,
      "sample_metadata_maaslin2.tsv"
    )
  )

  settings <- tibble(
    parameter = c(
      "Analysis",
      "Samples",
      "Viral families",
      "Normalization",
      "Transformation",
      "Analysis method",
      "Reference group",
      "Fixed effects",
      "Multiple-testing correction",
      "Significance threshold"
    ),
    value = c(
      analysis_name,
      nrow(feature_table),
      ncol(feature_table),
      "TSS",
      "LOG",
      "LM",
      "hibernation",
      model_formula,
      "BH",
      "q < 0.05"
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
  # MaAsLin2 model
  # -----------------------------------------------------------------------

  set.seed(12345)

  fit <- Maaslin2::Maaslin2(
    input_data = feature_table,
    input_metadata = metadata_table,
    output = model_dir,

    fixed_effects = fixed_effects,
    random_effects = NULL,

    reference = c(
      "group,hibernation"
    ),

    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",

    # Filtering was already performed before MaAsLin3.
    min_abundance = 0,
    min_prevalence = 0,
    min_variance = 0,

    correction = "BH",
    max_significance = 0.05,

    standardize = TRUE,

    # Custom visualization will be created later.
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    save_scatter = FALSE,
    save_models = FALSE,

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
      !is.na(qval),
      qval < 0.05
    ) |>
    arrange(
      comparison,
      qval
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

  model_summaries[[analysis_name]] <- tibble(
    analysis = analysis_name,
    samples = nrow(feature_table),
    ABH_samples = sum(metadata$group == "active"),
    H_samples = sum(metadata$group == "hibernation"),
    AAH_samples = sum(
      metadata$group == "post_hibernation"
    ),
    viral_families = ncol(feature_table),
    significant_group_associations =
      nrow(significant_group_results),
    significant_families =
      n_distinct(significant_group_results$feature),
    formula = model_formula
  )

  all_model_results[[analysis_name]] <- all_results
  all_group_results[[analysis_name]] <- group_results
  all_significant_results[[analysis_name]] <-
    significant_group_results
}

# -------------------------------------------------------------------------
# Combine all models
# -------------------------------------------------------------------------

model_summary <- bind_rows(
  model_summaries
)

combined_all_results <- bind_rows(
  all_model_results
)

combined_group_results <- bind_rows(
  all_group_results
)

combined_significant_results <- bind_rows(
  all_significant_results
)

write_tsv(
  model_summary,
  file.path(
    analysis_root,
    "model_summary.tsv"
  )
)

write_tsv(
  combined_all_results,
  file.path(
    analysis_root,
    "all_models_all_results.tsv"
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

cat(
  "\nMaAsLin2 method-sensitivity analysis completed.\n",
  "Results: ",
  analysis_root,
  "\n",
  sep = ""
)

print(model_summary)
