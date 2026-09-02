# Visualize BLASTn HSPs relative to the Previous_P5 reference.

library(ggplot2)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript 11_plot_densovirus_structure.R",
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

plot_data <- blast_hits %>%
  filter(
    query_id != "Previous_P5",
    alignment_length >= 50
  ) %>%
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
        status,
        length_nt
      ),
    by = c(
      "query_id" = "analysis_id"
    )
  )

orientation_summary <- plot_data %>%
  group_by(
    query_id,
    dataset,
    group,
    status,
    length_nt
  ) %>%
  summarise(
    hsp_count = n(),
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
  )

write_tsv(
  orientation_summary,
  file.path(
    output_directory,
    "densovirus_orientation_summary.tsv"
  )
)

make_structure_plot <- function(
  current_data,
  plot_title,
  output_file,
  number_columns,
  plot_width,
  plot_height
) {

  structure_plot <- ggplot(
    current_data
  ) +
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
      linewidth = 1.2,
      alpha = 0.9
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
      title = plot_title,
      subtitle = paste(
        "Each line represents one BLASTn HSP;",
        "Previous_P5 coordinates are shown on the y-axis"
      ),
      x = "Position in query sequence (nt)",
      y = "Position in Previous_P5 (nt)",
      colour = "Orientation"
    ) +
    theme_classic(
      base_size = 13
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      strip.text = element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    )

  ggsave(
    filename = output_file,
    plot = structure_plot,
    width = plot_width,
    height = plot_height,
    dpi = 300
  )
}

current_study_data <- plot_data %>%
  filter(
    dataset == "current_study"
  )

previous_study_data <- plot_data %>%
  filter(
    dataset == "previous_study"
  )

make_structure_plot(
  current_data = current_study_data,
  plot_title = paste(
    "Current near-complete densovirus sequences",
    "compared with Previous_P5"
  ),
  output_file = file.path(
    output_directory,
    "current_densoviruses_vs_previous_p5.png"
  ),
  number_columns = 2,
  plot_width = 12,
  plot_height = 14
)

make_structure_plot(
  current_data = previous_study_data,
  plot_title = paste(
    "Previous-study densovirus genomes",
    "compared with Previous_P5"
  ),
  output_file = file.path(
    output_directory,
    "previous_densoviruses_vs_previous_p5.png"
  ),
  number_columns = 3,
  plot_width = 15,
  plot_height = 12
)

message("Densovirus structure visualization completed.")
message("Results: ", output_directory)
