# Quality control of PICRUSt2 functional predictions.

library(dplyr)
library(readr)
library(ggplot2)

source("workflows/picrust2/r/local_paths.R")

# Input files.
fasta_file <- file.path(
  picrust2_input_dir,
  "input_asv_sequences.fasta"
)

asv_nsti_file <- file.path(
  picrust2_input_dir,
  "asv_nsti.tsv.gz"
)

weighted_nsti_file <- file.path(
  picrust2_input_dir,
  "sample_weighted_nsti.tsv.gz"
)

metadata_file <- file.path(
  picrust2_input_dir,
  "sample_metadata.tsv"
)

asv_counts_file <- file.path(
  picrust2_input_dir,
  "asv_counts_unrarefied.tsv"
)

# Output directory.
qc_output_dir <- file.path(
  picrust2_results_dir,
  "qc"
)

dir.create(
  qc_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Experimental groups.
group_levels <- c(
  "active",
  "hibernation",
  "post_hibernation"
)

group_labels <- c(
  active = "ABH",
  hibernation = "H",
  post_hibernation = "AAH"
)

group_colors <- c(
  active = "#0000FF",
  hibernation = "#FF0000",
  post_hibernation = "gray"
)

# Read identifiers of all input ASVs from FASTA.
fasta_lines <- readLines(fasta_file)

input_asv_ids <- fasta_lines[
  grepl("^>", fasta_lines)
] |>
  sub(
    pattern = "^>",
    replacement = ""
  ) |>
  sub(
    pattern = "\\s.*$",
    replacement = ""
  )

# Read the ASV-level NSTI table.
asv_nsti <- read_tsv(
  asv_nsti_file,
  show_col_types = FALSE
)

# Read sample-level weighted NSTI.
weighted_nsti <- read_tsv(
  weighted_nsti_file,
  show_col_types = FALSE
) |>
  rename(sample_id = sample)

# Read biological sample metadata.
metadata <- read_tsv(
  metadata_file,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  rename(sample_id = 1) |>
  mutate(
    group = factor(
      group,
      levels = group_levels
    )
  )

# Read the original unrarefied ASV counts.
asv_counts <- read_tsv(
  asv_counts_file,
  skip = 1,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  rename(feature_id = 1)

# Identify ASVs that were not placed by PICRUSt2.
placed_asv_ids <- asv_nsti$sequence

unplaced_asv_ids <- setdiff(
  input_asv_ids,
  placed_asv_ids
)

# Summarize the abundance of unplaced ASVs.
unplaced_counts <- asv_counts |>
  filter(feature_id %in% unplaced_asv_ids)

if (nrow(unplaced_counts) > 0) {

  unplaced_matrix <- unplaced_counts |>
    select(-feature_id) |>
    as.matrix()

  unplaced_summary <- tibble(
    feature_id = unplaced_counts$feature_id,
    total_reads = rowSums(unplaced_matrix),
    positive_samples = rowSums(unplaced_matrix > 0),
    maximum_sample_count = apply(
      unplaced_matrix,
      1,
      max
    )
  )

} else {

  unplaced_summary <- tibble(
    feature_id = character(),
    total_reads = numeric(),
    positive_samples = numeric(),
    maximum_sample_count = numeric()
  )
}

write_tsv(
  unplaced_summary,
  file.path(
    qc_output_dir,
    "unplaced_asvs.tsv"
  )
)

# Join weighted NSTI with metadata.
sample_qc <- metadata |>
  inner_join(
    weighted_nsti,
    by = "sample_id"
  )

if (nrow(sample_qc) != nrow(metadata)) {
  stop(
    "Not all metadata samples have weighted NSTI values."
  )
}

write_tsv(
  sample_qc |>
    mutate(group = as.character(group)),
  file.path(
    qc_output_dir,
    "sample_weighted_nsti.tsv"
  )
)

# Summarize weighted NSTI within each experimental group.
weighted_nsti_by_group <- sample_qc |>
  group_by(group) |>
  summarise(
    samples = n(),
    mean_weighted_nsti = mean(weighted_NSTI),
    median_weighted_nsti = median(weighted_NSTI),
    minimum_weighted_nsti = min(weighted_NSTI),
    maximum_weighted_nsti = max(weighted_NSTI),
    .groups = "drop"
  ) |>
  mutate(group = as.character(group))

write_tsv(
  weighted_nsti_by_group,
  file.path(
    qc_output_dir,
    "weighted_nsti_by_group.tsv"
  )
)

# Check whether weighted NSTI differs among groups.
weighted_nsti_test <- kruskal.test(
  weighted_NSTI ~ group,
  data = sample_qc
)

weighted_nsti_test_table <- tibble(
  test = "Kruskal-Wallis",
  statistic = unname(weighted_nsti_test$statistic),
  degrees_of_freedom = unname(weighted_nsti_test$parameter),
  p_value = weighted_nsti_test$p.value
)

write_tsv(
  weighted_nsti_test_table,
  file.path(
    qc_output_dir,
    "weighted_nsti_kruskal_wallis.tsv"
  )
)

# Create the overall PICRUSt2 QC summary.
qc_summary <- tibble(
  metric = c(
    "Input ASVs",
    "Placed ASVs",
    "Unplaced ASVs",
    "Placed ASVs (%)",
    "Bacterial placements",
    "Archaeal placements",
    "Mean ASV NSTI",
    "Median ASV NSTI",
    "Maximum ASV NSTI",
    "ASVs with NSTI > 0.15",
    "ASVs with NSTI > 0.30",
    "ASVs with NSTI > 1",
    "ASVs with NSTI > 2",
    "Samples",
    "Mean weighted NSTI",
    "Median weighted NSTI",
    "Minimum weighted NSTI",
    "Maximum weighted NSTI"
  ),
  value = c(
    length(input_asv_ids),
    length(placed_asv_ids),
    length(unplaced_asv_ids),
    100 * length(placed_asv_ids) /
      length(input_asv_ids),
    sum(asv_nsti$best_domain == "bac"),
    sum(asv_nsti$best_domain == "arc"),
    mean(asv_nsti$metadata_NSTI),
    median(asv_nsti$metadata_NSTI),
    max(asv_nsti$metadata_NSTI),
    sum(asv_nsti$metadata_NSTI > 0.15),
    sum(asv_nsti$metadata_NSTI > 0.30),
    sum(asv_nsti$metadata_NSTI > 1),
    sum(asv_nsti$metadata_NSTI > 2),
    nrow(sample_qc),
    mean(sample_qc$weighted_NSTI),
    median(sample_qc$weighted_NSTI),
    min(sample_qc$weighted_NSTI),
    max(sample_qc$weighted_NSTI)
  )
)

write_tsv(
  qc_summary,
  file.path(
    qc_output_dir,
    "picrust2_qc_summary.tsv"
  )
)

# Format the QC test result for the figure.
p_label <- if (
  weighted_nsti_test$p.value < 0.001
) {
  "K–W, p < 0.001"
} else {
  paste0(
    "K–W, p = ",
    formatC(
      weighted_nsti_test$p.value,
      digits = 3,
      format = "f"
    )
  )
}

# Plot weighted NSTI by experimental group.
weighted_nsti_plot <- ggplot(
  sample_qc,
  aes(
    x = group,
    y = weighted_NSTI,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.12,
    shape = 21,
    size = 2.7,
    color = "black"
  ) +
  annotate(
    geom = "text",
    x = 1,
    y = Inf,
    label = p_label,
    hjust = 0,
    vjust = 1.5,
    size = 5
  ) +
  scale_fill_manual(
    values = group_colors,
    guide = "none"
  ) +
  scale_x_discrete(
    labels = group_labels,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Weighted NSTI"
  ) +
  theme_classic(
    base_size = 16
  ) +
  theme(
    axis.text.x = element_text(
      size = 15
    ),
    axis.title.y = element_text(
      size = 17
    )
  )

ggsave(
  filename = file.path(
    qc_output_dir,
    "weighted_nsti_by_group.png"
  ),
  plot = weighted_nsti_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("")
message("PICRUSt2 QC completed.")
message("Results: ", qc_output_dir)
message("")
print(qc_summary)
message("")
print(weighted_nsti_by_group)
message("")
print(weighted_nsti_test_table)
