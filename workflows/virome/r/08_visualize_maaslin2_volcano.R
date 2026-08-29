# Visualize the exploratory MaAsLin2 analysis as a volcano plot.

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggtext)

source("workflows/virome/r/local_paths.R")

# Input and output paths.
maaslin2_dir <- file.path(
  virome_results_dir,
  "maaslin2",
  "viral_families_method_sensitivity"
)

input_file <- file.path(
  maaslin2_dir,
  "primary_all_samples",
  "group_results.tsv"
)

figure_dir <- file.path(maaslin2_dir, "figures")

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Visualization thresholds.
q_threshold <- 0.05
beta_threshold <- 2

# Read MaAsLin2 results.
results <- read_tsv(
  input_file,
  show_col_types = FALSE
) %>%
  filter(metadata == "group") %>%
  mutate(
    qval = pmax(qval, 1e-300),

    significant = (
      qval < q_threshold &
      abs(coef) >= beta_threshold
    ),

    comparison_class = case_when(
      significant & value == "active" ~ "ABH vs H",
      significant & value == "post_hibernation" ~ "AAH vs H",
      TRUE ~ "Not significant"
    ),

    comparison_class = factor(
      comparison_class,
      levels = c(
        "ABH vs H",
        "AAH vs H",
        "Not significant"
      )
    ),

    taxon_label = if_else(
      significant,
      feature,
      NA_character_
    ),

    minus_log10_q = -log10(qval)
  )

# Save the data used for the figure.
write_tsv(
  results,
  file.path(
    figure_dir,
    "maaslin2_viral_family_volcano_primary_data.tsv"
  )
)

# Colours matching the previous virome figure.
volcano_colours <- c(
  "ABH vs H" = "red",
  "AAH vs H" = "blue",
  "Not significant" = "#BDBDBD"
)

volcano_plot <- ggplot(
  results,
  aes(
    x = coef,
    y = minus_log10_q,
    colour = comparison_class
  )
) +
  geom_hline(
    yintercept = -log10(q_threshold),
    colour = "grey70",
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = c(-beta_threshold, beta_threshold),
    colour = "grey70",
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  geom_point(
    size = 2.8,
    alpha = 0.95
  ) +
  geom_text_repel(
    aes(label = taxon_label),
    fontface = "italic",
    size = 5.2,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE,
    seed = 123
  ) +
  scale_colour_manual(
    values = volcano_colours,
    breaks = c(
      "ABH vs H",
      "AAH vs H",
      "Not significant"
    ),
    labels = c(
      "<b>ABH vs H</b>:<br><i>q</i> < 0.05 and |β| ≥ 2",
      "<b>AAH vs H</b>:<br><i>q</i> < 0.05 and |β| ≥ 2",
      "<b>Not significant</b>"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    labels = function(x) {
      gsub("-", "\u2212", format(x, trim = TRUE))
    },
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.12))
  ) +
  labs(
    x = "β coefficient",
    y = "−log<sub>10</sub>(<i>q</i>)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),

    axis.title.x = element_text(
      size = 19,
      margin = margin(t = 10)
    ),

    axis.title.y = element_markdown(
      size = 19,
      margin = margin(r = 10)
    ),

    axis.text = element_text(size = 14),

    legend.position = "right",
    legend.text = element_markdown(
      size = 15,
      lineheight = 1.1
    ),

    legend.key.height = grid::unit(0.75, "cm"),

    plot.margin = margin(
      t = 15,
      r = 20,
      b = 15,
      l = 15
    )
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

# Save PNG and vector PDF versions.
ggsave(
  filename = file.path(
    figure_dir,
    "maaslin2_viral_family_volcano_primary.png"
  ),
  plot = volcano_plot,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_dir,
    "maaslin2_viral_family_volcano_primary.pdf"
  ),
  plot = volcano_plot,
  width = 12,
  height = 7,
  device = cairo_pdf
)

cat(
  "MaAsLin2 volcano plot completed.\n",
  "Results:", figure_dir, "\n",
  "Highlighted associations:",
  sum(results$significant),
  "\n"
)
