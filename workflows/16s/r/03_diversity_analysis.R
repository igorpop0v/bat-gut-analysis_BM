# Visualize the selected alpha- and beta-diversity metrics.
# Statistical results are read from existing QIIME 2 visualizations.

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(patchwork)

source("workflows/16s/r/local_paths.R")

qiime2_dir <- dirname(input_dir)

output_dir <- file.path(
  results_dir,
  "diversity"
)

dir.create(
  output_dir,
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

status_labels <- c(
  active = "Before hibernation",
  hibernation = "Hibernation",
  post_hibernation = "After hibernation"
)

# Find a file stored inside a QIIME 2 visualization.
qzv_member <- function(qzv_file, ending) {

  members <- unzip(
    qzv_file,
    list = TRUE
  )$Name

  member <- members[
    endsWith(members, ending)
  ]

  if (length(member) != 1) {
    stop(
      "Could not find ",
      ending,
      " in ",
      qzv_file
    )
  }

  member
}


# Read Kruskal-Wallis statistics from one QIIME 2 visualization.
read_alpha_statistics <- function(qzv_file, metric) {

  json_file <- qzv_member(
    qzv_file,
    "/data/column-group.jsonp"
  )

  json_text <- paste(
    readLines(
      unz(qzv_file, json_file),
      warn = FALSE
    ),
    collapse = ""
  )

  overall_values <- str_match(
    json_text,
    '"H": ([^,]+), "p": ([^}]+)'
  )

  pairwise_file <- qzv_member(
    qzv_file,
    "/data/kruskal-wallis-pairwise-group.csv"
  )

  pairwise <- read_csv(
    unz(qzv_file, pairwise_file),
    show_col_types = FALSE
  ) |>
    rename(
      group_1 = `Group 1`,
      group_2 = `Group 2`,
      p_value = `p-value`,
      q_value = `q-value`
    ) |>
    mutate(
      metric = metric,
      group_1 = str_remove(
        group_1,
        " \\(n=.*$"
      ),
      group_2 = str_remove(
        group_2,
        " \\(n=.*$"
      )
    ) |>
    select(
      metric,
      group_1,
      group_2,
      H,
      p_value,
      q_value
    )

  list(
    overall = tibble(
      metric = metric,
      H = as.numeric(
        overall_values[1, 2]
      ),
      p_value = as.numeric(
        overall_values[1, 3]
      )
    ),
    pairwise = pairwise
  )
}


# Read the overall PERMANOVA or PERMDISP result.
read_beta_overall <- function(qzv_file) {

  html_file <- qzv_member(
    qzv_file,
    "/data/index.html"
  )

  html <- paste(
    readLines(
      unz(qzv_file, html_file),
      warn = FALSE
    ),
    collapse = " "
  )

  values <- str_match_all(
    html,
    "<td>(.*?)</td>"
  )[[1]][, 2] |>
    str_squish()

  tibble(
    method = values[1],
    statistic_name = values[2],
    statistic = as.numeric(values[5]),
    p_value = as.numeric(values[6]),
    permutations = as.integer(values[7])
  )
}


# Read the first two PCoA axes.
read_pcoa <- function(pcoa_file) {

  lines <- readLines(
    pcoa_file,
    warn = FALSE
  )

  proportion_row <- grep(
    "^Proportion explained\\t",
    lines
  )

  proportion <- as.numeric(
    str_split(
      lines[proportion_row + 1],
      "\\t"
    )[[1]]
  )

  site_row <- grep(
    "^Site\\t",
    lines
  )

  number_of_sites <- as.integer(
    str_split(
      lines[site_row],
      "\\t"
    )[[1]][2]
  )

  site_lines <- lines[
    site_row + seq_len(number_of_sites)
  ]

  coordinates <- read_tsv(
    I(paste(site_lines, collapse = "\n")),
    col_names = FALSE,
    show_col_types = FALSE
  ) |>
    transmute(
      sample_id = X1,
      PCoA1 = X2,
      PCoA2 = X3
    )

  list(
    coordinates = coordinates,
    proportion = proportion
  )
}


# Read metadata.
metadata <- read_tsv(
  file.path(
    input_dir,
    "sample_metadata.tsv"
  ),
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


# Read alpha-diversity values.
alpha <- read_tsv(
  file.path(
    input_dir,
    "alpha_diversity_metrics_rarefied_18000.tsv"
  ),
  show_col_types = FALSE
) |>
  left_join(
    metadata,
    by = "sample_id"
  ) |>
  pivot_longer(
    c(
      shannon,
      simpson,
      pielou_evenness
    ),
    names_to = "metric",
    values_to = "value"
  )


# QIIME 2 alpha-diversity visualizations.
alpha_qzv <- c(
  shannon = "shannon_group_significance.qzv",
  simpson = "simpson_group_significance.qzv",
  pielou_evenness = "evenness_group_significance.qzv"
)

alpha_statistics <- lapply(
  names(alpha_qzv),
  function(metric) {

    read_alpha_statistics(
      file.path(
        qiime2_dir,
        "diversity",
        "alpha_significance",
        alpha_qzv[[metric]]
      ),
      metric
    )
  }
)

alpha_overall <- bind_rows(
  lapply(
    alpha_statistics,
    `[[`,
    "overall"
  )
)

alpha_pairwise <- bind_rows(
  lapply(
    alpha_statistics,
    `[[`,
    "pairwise"
  )
)


# QIIME 2 beta-diversity statistics.
permanova_overall <- read_beta_overall(
  file.path(
    qiime2_dir,
    "diversity",
    "bray_curtis_significance",
    "bray_curtis_permanova.qzv"
  )
)

permdisp_overall <- read_beta_overall(
  file.path(
    qiime2_dir,
    "diversity",
    "bray_curtis_significance",
    "bray_curtis_permdisp.qzv"
  )
)

permanova_pairwise <- read_csv(
  file.path(
    input_dir,
    "bray_curtis_permanova_pairwise.csv"
  ),
  show_col_types = FALSE
)

permdisp_pairwise <- read_csv(
  file.path(
    input_dir,
    "bray_curtis_permdisp_pairwise.csv"
  ),
  show_col_types = FALSE
)


# Save all statistical results as ordinary tables.
write_tsv(
  alpha_overall,
  file.path(
    output_dir,
    "alpha_kruskal_wallis_overall.tsv"
  )
)

write_tsv(
  alpha_pairwise,
  file.path(
    output_dir,
    "alpha_kruskal_wallis_pairwise.tsv"
  )
)

write_tsv(
  permanova_overall,
  file.path(
    output_dir,
    "bray_curtis_permanova_overall.tsv"
  )
)

write_tsv(
  permdisp_overall,
  file.path(
    output_dir,
    "bray_curtis_permdisp_overall.tsv"
  )
)

write_csv(
  permanova_pairwise,
  file.path(
    output_dir,
    "bray_curtis_permanova_pairwise.csv"
  )
)

write_csv(
  permdisp_pairwise,
  file.path(
    output_dir,
    "bray_curtis_permdisp_pairwise.csv"
  )
)


# Format p- and q-values for the figure.
format_p <- function(p) {

  if (p < 0.001) {
    "p < 0.001"
  } else {
    paste0(
      "p = ",
      formatC(
        p,
        digits = 3,
        format = "f"
      )
    )
  }
}

format_q <- function(q) {

  if (q < 0.001) {
    "q < 0.001"
  } else {
    paste0(
      "q = ",
      formatC(
        q,
        digits = 3,
        format = "f"
      )
    )
  }
}


# Create one alpha-diversity panel.
make_alpha_plot <- function(metric_name, y_title) {
  
  plot_data <- alpha |>
    filter(metric == metric_name)
  
  overall <- alpha_overall |>
    filter(metric == metric_name)
  
  comparisons <- alpha_pairwise |>
    filter(
      metric == metric_name,
      q_value < 0.05
    ) |>
    mutate(
      x1 = match(
        group_1,
        group_levels
      ),
      x2 = match(
        group_2,
        group_levels
      )
    )
  
  y_min <- min(plot_data$value)
  y_max <- max(plot_data$value)
  y_range <- y_max - y_min
  
  # Put every pairwise comparison at a different height.
  comparisons <- comparisons |>
    arrange(x1, x2) |>
    mutate(
      bracket_y =
        y_max +
        (
          0.08 +
            0.14 * (row_number() - 1)
        ) * y_range,
      
      text_y =
        bracket_y +
        0.055 * y_range
    )
  
  overall_y <- y_max + 0.42 * y_range
  upper_limit <- y_max + 0.48 * y_range
  
  ggplot(
    plot_data,
    aes(
      group,
      value,
      fill = group
    )
  ) +
    
    # Median line.
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.55,
      linewidth = 0.8,
      color = "black",
      show.legend = FALSE
    ) +
    
    # Individual samples.
    geom_point(
      shape = 21,
      color = "black",
      size = 2.5,
      stroke = 0.55,
      position = position_jitter(
        width = 0.22,
        height = 0,
        seed = 10000
      ),
      show.legend = FALSE
    ) +
    
    # Pairwise-comparison lines.
    geom_segment(
      data = comparisons,
      aes(
        x = x1,
        xend = x2,
        y = bracket_y,
        yend = bracket_y
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.6
    ) +
    
    # Pairwise q-values above the lines.
    geom_text(
      data = comparisons,
      aes(
        x = (x1 + x2) / 2,
        y = text_y,
        label = vapply(
          q_value,
          format_q,
          character(1)
        )
      ),
      inherit.aes = FALSE,
      size = 4
    ) +
    
    # Overall Kruskal-Wallis result.
    annotate(
      "text",
      x = 1,
      y = overall_y,
      label = paste0(
        "K-W, ",
        format_p(overall$p_value)
      ),
      hjust = 0,
      size = 4.2
    ) +
    
    scale_fill_manual(
      values = group_colors,
      guide = "none"
    ) +
    
    scale_y_continuous(
      limits = c(
        y_min,
        upper_limit
      )
    ) +
    
    labs(
      x = NULL,
      y = y_title
    ) +
    
    theme_classic(
      base_size = 16
    ) +
    
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

plot_shannon <- make_alpha_plot(
  "shannon",
  "Shannon index"
)

plot_simpson <- make_alpha_plot(
  "simpson",
  "Simpson index"
)

plot_pielou <- make_alpha_plot(
  "pielou_evenness",
  "Pielou evenness"
)


# Read the Bray-Curtis PCoA coordinates.
pcoa <- read_pcoa(
  file.path(
    input_dir,
    "bray_curtis_pcoa_ordination_rarefied_18000.txt"
  )
)

pcoa_data <- pcoa$coordinates |>
  left_join(
    metadata,
    by = "sample_id"
  )


# Statistical annotation for the PCoA panel.
beta_label <- paste0(
  "PERMANOVA: pseudo-F = ",
  formatC(
    permanova_overall$statistic,
    digits = 2,
    format = "f"
  ),
  ", ",
  format_p(permanova_overall$p_value),
  "\nPairwise PERMANOVA: all q ≤ ",
  formatC(
    max(permanova_pairwise$`q-value`),
    digits = 4,
    format = "f"
  ),
  "\nPERMDISP: F = ",
  formatC(
    permdisp_overall$statistic,
    digits = 2,
    format = "f"
  ),
  ", ",
  format_p(permdisp_overall$p_value)
)


# Bray-Curtis PCoA.
plot_pcoa <- ggplot(
  pcoa_data,
  aes(
    PCoA1,
    PCoA2,
    color = group
  )
) +
  
  # Samples are shown as simple colored points.
  geom_point(
    size = 2.5
  ) +
  
  # Ellipses are contours only, without any fill.
  stat_ellipse(
    geom = "path",
    level = 0.95,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  # Statistics are placed in the additional space above the plot.
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = beta_label,
    hjust = -0.05,
    vjust = 1.05,
    size = 4.2,
    fill = "white",
    linewidth = 0
  ) +
  
  scale_color_manual(
    values = group_colors,
    labels = group_labels,
    name = "Group"
  ) +
  
  # Display proper minus signs.
  scale_x_continuous(
    labels = function(x) {
      sub(
        "-",
        "−",
        x,
        fixed = TRUE
      )
    },
    expand = expansion(
      mult = c(0.05, 0.05)
    )
  ) +
  
  # The larger upper expansion creates space for the statistics.
  scale_y_continuous(
    labels = function(x) {
      sub(
        "-",
        "−",
        x,
        fixed = TRUE
      )
    },
    expand = expansion(
      mult = c(0.05, 0.38)
    )
  ) +
  
  labs(
    x = paste0(
      "PCoA axis 1 (",
      round(
        100 * pcoa$proportion[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PCoA axis 2 (",
      round(
        100 * pcoa$proportion[2],
        1
      ),
      "%)"
    )
  ) +
  
  theme_classic(
    base_size = 16
  ) +
  
  theme(
    legend.title = element_text(
      size = 16,
      face = "plain"
    ),
    legend.text = element_text(
      size = 16
    ),
    legend.key = element_blank()
  ) +
  
  # The legend contains points only.
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 16,
        size = 3,
        linewidth = 0
      )
    )
  )


# Combine all panels and collect the legend.
final_figure <- (
  plot_shannon |
    plot_simpson |
    plot_pielou
) /
  plot_pcoa +
  
  plot_layout(
    heights = c(1, 1.15),
    guides = "collect"
  ) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    legend.position = "right"
  )


# Save PNG and vector PDF versions.
ggsave(
  file.path(
    output_dir,
    "alpha_beta_diversity.png"
  ),
  final_figure,
  width = 15,
  height = 10,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "alpha_beta_diversity.pdf"
  ),
  final_figure,
  width = 15,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

message(
  "Diversity analysis complete. Results: ",
  output_dir
)
