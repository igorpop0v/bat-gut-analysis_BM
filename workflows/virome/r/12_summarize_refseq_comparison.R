# Summarize BLASTn comparison with the public RefSeq densovirus genome.

library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript 12_summarize_refseq_comparison.R",
      "<blast_hsps.tsv> <manifest.tsv> <output_directory>"
    )
  )
}

blast_file <- args[1]
manifest_file <- args[2]
output_directory <- args[3]

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

column_names <- c(
  "query_id",
  "reference_id",
  "identity_pct",
  "alignment_length",
  "query_start",
  "query_end",
  "query_length",
  "reference_start",
  "reference_end",
  "reference_length",
  "e_value",
  "bitscore"
)

blast_hits <- read_tsv(
  blast_file,
  col_names = column_names,
  show_col_types = FALSE
)

manifest <- read_tsv(
  manifest_file,
  show_col_types = FALSE
)

count_covered_bases <- function(starts, ends) {
  positions <- Map(
    function(start, end) {
      seq.int(min(start, end), max(start, end))
    },
    starts,
    ends
  )

  length(unique(unlist(positions, use.names = FALSE)))
}

comparison_data <- blast_hits %>%
  filter(alignment_length >= 50) %>%
  mutate(
    hsp_orientation = if_else(
      reference_start <= reference_end,
      "Forward",
      "Reverse"
    )
  ) %>%
  left_join(
    manifest %>%
      select(
        analysis_id,
        dataset,
        group,
        status
      ),
    by = c("query_id" = "analysis_id")
  )

comparison_summary <- comparison_data %>%
  group_by(
    query_id,
    dataset,
    group,
    status
  ) %>%
  summarise(
    query_length = first(query_length),
    reference_length = first(reference_length),
    hsp_count = n(),
    best_hsp_identity_pct = max(identity_pct),
    weighted_hsp_identity_pct = weighted.mean(
      identity_pct,
      alignment_length
    ),
    query_covered_bases = count_covered_bases(
      query_start,
      query_end
    ),
    reference_covered_bases = count_covered_bases(
      reference_start,
      reference_end
    ),
    forward_bitscore = sum(
      bitscore[hsp_orientation == "Forward"]
    ),
    reverse_bitscore = sum(
      bitscore[hsp_orientation == "Reverse"]
    ),
    dominant_orientation = if_else(
      forward_bitscore >= reverse_bitscore,
      "Forward",
      "Reverse"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    query_coverage_pct =
      100 * query_covered_bases / query_length,
    reference_coverage_pct =
      100 * reference_covered_bases / reference_length
  ) %>%
  arrange(dataset, query_id)

write_tsv(
  comparison_summary,
  file.path(
    output_directory,
    "densovirus_vs_refseq_summary.tsv"
  )
)

orientation_summary <- comparison_summary %>%
  select(
    query_id,
    dataset,
    group,
    status,
    query_length,
    forward_bitscore,
    reverse_bitscore,
    dominant_orientation
  )

write_tsv(
  orientation_summary,
  file.path(
    output_directory,
    "densovirus_refseq_orientation_summary.tsv"
  )
)

make_structure_plot <- function(
  data,
  title,
  filename,
  number_columns,
  width,
  height
) {
  plot <- ggplot(data) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linewidth = 0.5,
      linetype = "dashed",
      colour = "grey70"
    ) +
    geom_segment(
      aes(
        x = query_start,
        xend = query_end,
        y = reference_start,
        yend = reference_end,
        colour = hsp_orientation
      ),
      linewidth = 1.1
    ) +
    facet_wrap(
      ~ query_id,
      scales = "free_x",
      ncol = number_columns
    ) +
    scale_colour_manual(
      values = c(
        "Forward" = "#0072B2",
        "Reverse" = "#D55E00"
      )
    ) +
    labs(
      title = title,
      x = "Position in query sequence (nt)",
      y = "Position in RefSeq NC_031450.1 (nt)",
      colour = "Orientation"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(
    filename = file.path(output_directory, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

make_structure_plot(
  comparison_data %>%
    filter(dataset == "current_study"),
  "Current densovirus sequences compared with RefSeq NC_031450.1",
  "current_densoviruses_vs_refseq.png",
  2,
  12,
  14
)

make_structure_plot(
  comparison_data %>%
    filter(dataset == "previous_study"),
  "Previous-study densovirus sequences compared with RefSeq NC_031450.1",
  "previous_densoviruses_vs_refseq.png",
  3,
  15,
  12
)

message("RefSeq comparison completed.")
message("Results: ", output_directory)
