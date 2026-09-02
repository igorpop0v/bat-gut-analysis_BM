# Summarize and visualize read depth across densovirus candidates.

library(ggplot2)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript 10_plot_densovirus_read_depth.R",
      "<mapping_summary.tsv> <depth_directory> <output_directory>"
    )
  )
}

mapping_summary_file <- args[1]
depth_directory <- args[2]
output_directory <- args[3]

figure_directory <- file.path(
  output_directory,
  "figures"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

mapping_summary <- read_tsv(
  mapping_summary_file,
  show_col_types = FALSE
)

safe_filename <- function(candidate_id) {

  gsub(
    "[^A-Za-z0-9._-]",
    "_",
    candidate_id
  )
}

summary_list <- vector(
  "list",
  nrow(mapping_summary)
)

plot_list <- vector(
  "list",
  nrow(mapping_summary)
)

for (i in seq_len(nrow(mapping_summary))) {

  candidate_id <- mapping_summary$candidate_id[i]
  sample_id <- mapping_summary$sample_id[i]

  depth_filename <- paste0(
    gsub("\\|", "__", candidate_id),
    ".depth.tsv"
  )

  depth_file <- file.path(
    depth_directory,
    depth_filename
  )

  depth_table <- read_tsv(
    depth_file,
    col_names = c(
      "contig",
      "position",
      "depth"
    ),
    show_col_types = FALSE
  )

  candidate_length <- max(
    depth_table$position
  )

  first_terminal <- depth_table$position <= 100

  second_terminal <- depth_table$position >
    candidate_length - 100

  terminal_positions <- first_terminal |
    second_terminal

  internal_positions <- depth_table$position > 100 &
    depth_table$position <= candidate_length - 100

  terminal_mean_depth <- mean(
    depth_table$depth[terminal_positions]
  )

  if (sum(internal_positions) > 0) {

    internal_mean_depth <- mean(
      depth_table$depth[internal_positions]
    )

    terminal_to_internal_ratio <-
      terminal_mean_depth /
      internal_mean_depth

  } else {

    internal_mean_depth <- NA_real_
    terminal_to_internal_ratio <- NA_real_
  }

  mean_depth <- mean(
    depth_table$depth
  )

  depth_cv <- ifelse(
    mean_depth > 0,
    sd(depth_table$depth) / mean_depth,
    NA_real_
  )

  summary_list[[i]] <- data.frame(
    candidate_id = candidate_id,
    sample_id = sample_id,
    candidate_length = candidate_length,
    minimum_depth = min(depth_table$depth),
    depth_percentile_5 = as.numeric(
      quantile(
        depth_table$depth,
        probs = 0.05,
        names = FALSE
      )
    ),
    median_depth = median(depth_table$depth),
    mean_depth = mean_depth,
    depth_percentile_95 = as.numeric(
      quantile(
        depth_table$depth,
        probs = 0.95,
        names = FALSE
      )
    ),
    maximum_depth = max(depth_table$depth),
    depth_cv = depth_cv,
    coverage_ge_1x_pct = 100 *
      mean(depth_table$depth >= 1),
    coverage_ge_10x_pct = 100 *
      mean(depth_table$depth >= 10),
    coverage_ge_100x_pct = 100 *
      mean(depth_table$depth >= 100),
    terminal_mean_depth = terminal_mean_depth,
    internal_mean_depth = internal_mean_depth,
    terminal_to_internal_ratio =
      terminal_to_internal_ratio,
    stringsAsFactors = FALSE
  )

  # Average depth within 50-nucleotide windows.
  binned_depth <- depth_table %>%
    mutate(
      position_bin = ceiling(position / 50)
    ) %>%
    group_by(position_bin) %>%
    summarise(
      position = mean(position),
      depth = mean(depth),
      .groups = "drop"
    )

  depth_table <- depth_table %>%
    mutate(
      log_depth = log10(depth + 1)
    )

  binned_depth <- binned_depth %>%
    mutate(
      log_depth = log10(depth + 1)
    )

  plot_subtitle <- paste0(
    sample_id,
    "; length = ",
    candidate_length,
    " nt; mean depth = ",
    format(
      mean_depth,
      digits = 4,
      scientific = FALSE,
      big.mark = ","
    ),
    "×"
  )

  coverage_plot <- ggplot(
    depth_table,
    aes(
      x = position,
      y = log_depth
    )
  ) +
    geom_line(
      linewidth = 0.25,
      colour = "grey70",
      alpha = 0.6
    ) +
    geom_line(
      data = binned_depth,
      aes(
        x = position,
        y = log_depth
      ),
      linewidth = 0.8,
      colour = "#0072B2"
    ) +
    labs(
      title = candidate_id,
      subtitle = plot_subtitle,
      x = "Position (nt)",
      y = expression(log[10] * "(read depth + 1)")
    ) +
    theme_classic(
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      plot.subtitle = element_text(
        size = 11
      ),
      axis.title = element_text(
        face = "bold"
      )
    )

  plot_list[[i]] <- coverage_plot

  ggsave(
    filename = file.path(
      figure_directory,
      paste0(
        safe_filename(candidate_id),
        "_depth.png"
      )
    ),
    plot = coverage_plot,
    width = 10,
    height = 5,
    dpi = 300
  )
}

depth_summary <- bind_rows(
  summary_list
)

numeric_columns <- names(depth_summary)[
  vapply(
    depth_summary,
    is.numeric,
    logical(1)
  )
]

depth_summary[numeric_columns] <- lapply(
  depth_summary[numeric_columns],
  round,
  digits = 3
)

write_tsv(
  depth_summary,
  file.path(
    output_directory,
    "densovirus_read_depth_summary.tsv"
  )
)

pdf(
  file.path(
    output_directory,
    "densovirus_read_depth_profiles.pdf"
  ),
  width = 10,
  height = 5
)

for (coverage_plot in plot_list) {
  print(coverage_plot)
}

dev.off()

message("Densovirus depth analysis completed.")
message("Results: ", output_directory)
