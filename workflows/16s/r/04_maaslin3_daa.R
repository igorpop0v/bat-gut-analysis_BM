# Differential abundance and prevalence analysis with MaAsLin3.
# Taxonomic levels: phylum, class, order, family, and genus.
# Reference group: hibernation.

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(maaslin3)

# Local input and output paths.
source("workflows/16s/r/local_paths.R")

counts_file <- file.path(
  input_dir,
  "asv_counts_unrarefied.tsv"
)

metadata_file <- file.path(
  input_dir,
  "sample_metadata.tsv"
)

taxonomy_file <- file.path(
  results_dir,
  "taxonomy",
  "asv_taxonomy_parsed_to_genus.tsv"
)

# All results from this analysis will be stored here.
analysis_root <- file.path(
  results_dir,
  "maaslin3",
  "all_ranks_primary"
)

prepared_dir <- file.path(
  analysis_root,
  "prepared_inputs"
)

# Stop before accidentally overwriting an existing analysis.
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

# Taxonomic ranks included in the analysis.
taxonomic_ranks <- c(
  "phylum",
  "class",
  "order",
  "family",
  "genus"
)

# A taxon must be detected in at least 10% of samples.
min_prevalence <- 0.10

# Read the unrarefied ASV count table.
asv_counts <- read_tsv(
  counts_file,
  skip = 1,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  rename(feature_id = 1)

# Read parsed GTDB taxonomy.
taxonomy <- read_tsv(
  taxonomy_file,
  show_col_types = FALSE
)

# Read sample metadata.
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

# Convert ASV counts to long format.
asv_counts_long <- asv_counts |>
  pivot_longer(
    cols = -feature_id,
    names_to = "sample_id",
    values_to = "count"
  ) |>
  # The positive control is excluded because it is absent
  # from the biological sample metadata.
  filter(sample_id %in% metadata$sample_id)

# Calculate the original sequencing depth.
library_sizes <- asv_counts_long |>
  group_by(sample_id) |>
  summarise(
    library_size = sum(count),
    .groups = "drop"
  )

# Prepare metadata used in every rank-specific model.
model_metadata <- metadata |>
  left_join(
    library_sizes,
    by = "sample_id"
  ) |>
  mutate(
    log2_library_size = log2(library_size)
  ) |>
  arrange(match(sample_id, metadata$sample_id))

metadata_table <- model_metadata |>
  select(
    sample_id,
    group,
    log2_library_size
  ) |>
  column_to_rownames("sample_id") |>
  as.data.frame()

# Save the metadata used by MaAsLin3.
write_tsv(
  model_metadata |>
    mutate(group = as.character(group)),
  file.path(
    prepared_dir,
    "sample_metadata_maaslin3.tsv"
  )
)

# Objects used to combine results from all ranks.
rank_summaries <- list()
all_results <- list()
significant_results <- list()

# Run the same analysis independently at every taxonomic rank.
for (rank in taxonomic_ranks) {

  message("")
  message("Preparing taxonomic rank: ", rank)

  # Aggregate ASV counts at the selected taxonomic rank.
  rank_counts <- asv_counts_long |>
    inner_join(
      taxonomy |>
        transmute(
          feature_id,
          taxon = .data[[rank]]
        ),
      by = "feature_id"
    ) |>
    group_by(sample_id, taxon) |>
    summarise(
      count = sum(count),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = taxon,
      values_from = count,
      values_fill = 0
    ) |>
    arrange(match(sample_id, model_metadata$sample_id))

  # Check sample order.
  if (!identical(
    rank_counts$sample_id,
    model_metadata$sample_id
  )) {
    stop(
      "Sample IDs do not match at rank: ",
      rank
    )
  }

  # Check that no reads were lost during taxonomic aggregation.
  observed_totals <- rowSums(
    as.data.frame(rank_counts[, -1])
  )

  expected_totals <- model_metadata$library_size

  if (!isTRUE(all.equal(
    observed_totals,
    expected_totals,
    check.attributes = FALSE
  ))) {
    stop(
      "Counts were lost during aggregation at rank: ",
      rank
    )
  }

  # Calculate prevalence before filtering.
  rank_prevalence <- rank_counts |>
    pivot_longer(
      cols = -sample_id,
      names_to = "taxon",
      values_to = "count"
    ) |>
    group_by(taxon) |>
    summarise(
      positive_samples = sum(count > 0),
      prevalence = mean(count > 0),
      total_count = sum(count),
      .groups = "drop"
    ) |>
    mutate(
      passes_prevalence_filter =
        prevalence >= min_prevalence
    ) |>
    arrange(
      desc(prevalence),
      desc(total_count)
    )

  # Save rank-specific counts before normalization.
  write_tsv(
    rank_counts,
    file.path(
      prepared_dir,
      paste0(
        rank,
        "_counts_unrarefied.tsv"
      )
    )
  )

  # Save the prevalence filtering audit.
  write_tsv(
    rank_prevalence,
    file.path(
      prepared_dir,
      paste0(
        rank,
        "_prevalence_before_filtering.tsv"
      )
    )
  )

  # Record how many taxa pass filtering.
  rank_summaries[[rank]] <- tibble(
    taxonomic_rank = rank,
    samples = nrow(rank_counts),
    taxa_before_filtering = ncol(rank_counts) - 1,
    taxa_after_filtering =
      sum(rank_prevalence$passes_prevalence_filter),
    taxa_removed =
      sum(!rank_prevalence$passes_prevalence_filter),
    minimum_prevalence = min_prevalence,
    minimum_positive_samples =
      ceiling(min_prevalence * nrow(rank_counts)),
    normalization = "TSS",
    transformation = "LOG"
  )

  # Convert sample IDs to row names for MaAsLin3.
  feature_table <- rank_counts |>
    column_to_rownames("sample_id") |>
    as.data.frame(check.names = FALSE)

  rank_output_dir <- file.path(
    analysis_root,
    rank
  )

  message("Running MaAsLin3 for rank: ", rank)
  message(
    "Taxa before filtering: ",
    ncol(feature_table)
  )
  message(
    "Taxa passing filtering: ",
    sum(rank_prevalence$passes_prevalence_filter)
  )

  # Run MaAsLin3.
  fit_out <- maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata_table,
    output = rank_output_dir,

    # ABH and AAH are compared with hibernation.
    formula = ~ group + log2_library_size,

    # Normalize the count table separately at each rank.
    normalization = "TSS",

    # Log2 transformation of positive relative abundances.
    transform = "LOG",

    # Retain taxa detected in at least 10% of samples.
    min_prevalence = min_prevalence,
    min_abundance = 0,
    zero_threshold = 0,

    # Account for compositionality in the abundance model.
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE,

    # Flag prevalence associations potentially induced
    # by stronger abundance associations.
    warn_prevalence = TRUE,

    # Stabilize logistic prevalence models.
    augment = TRUE,

    # Multiple-testing correction is performed separately
    # within each taxonomic rank.
    correction = "BH",
    max_significance = 0.05,

    # Standardize the continuous sequencing-depth covariate.
    standardize = TRUE,

    # Generate a standard summary plot for each rank.
    plot_summary_plot = TRUE,
    summary_plot_first_n = 25,

    # Individual association plots will be generated later
    # together with the custom visualization.
    plot_associations = FALSE,

    cores = 1
  )

  # Read the results and add the taxonomic rank.
  all_results[[rank]] <- read_tsv(
    file.path(
      rank_output_dir,
      "all_results.tsv"
    ),
    show_col_types = FALSE
  ) |>
    mutate(
      taxonomic_rank = rank,
      fdr_scope = "within_taxonomic_rank",
      .before = 1
    )

  significant_results[[rank]] <- read_tsv(
    file.path(
      rank_output_dir,
      "significant_results.tsv"
    ),
    show_col_types = FALSE
  ) |>
    mutate(
      taxonomic_rank = rank,
      fdr_scope = "within_taxonomic_rank",
      .before = 1
    )
}

# Combine the filtering summaries.
rank_filter_summary <- bind_rows(
  rank_summaries
)

write_tsv(
  rank_filter_summary,
  file.path(
    analysis_root,
    "rank_filter_summary.tsv"
  )
)

# Combine all rank-specific model results.
all_ranks_all <- bind_rows(
  all_results
)

all_ranks_significant <- bind_rows(
  significant_results
)

write_tsv(
  all_ranks_all,
  file.path(
    analysis_root,
    "all_ranks_all_results.tsv"
  )
)

write_tsv(
  all_ranks_all |>
    filter(metadata == "group"),
  file.path(
    analysis_root,
    "all_ranks_group_results.tsv"
  )
)

write_tsv(
  all_ranks_significant,
  file.path(
    analysis_root,
    "all_ranks_significant_results.tsv"
  )
)

write_tsv(
  all_ranks_significant |>
    filter(metadata == "group"),
  file.path(
    analysis_root,
    "all_ranks_group_significant_results.tsv"
  )
)

# Record package and R versions.
writeLines(
  capture.output(sessionInfo()),
  file.path(
    analysis_root,
    "session_info.txt"
  )
)

message("")
message("All taxonomic ranks completed successfully.")
message("Results: ", analysis_root)
