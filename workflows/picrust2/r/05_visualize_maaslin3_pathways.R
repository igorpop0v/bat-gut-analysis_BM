# Custom visualization of MaAsLin3 associations
# for predicted MetaCyc pathways.

library(dplyr)
library(readr)
library(ggplot2)
library(ggnewscale)
library(stringr)

source("workflows/picrust2/r/local_paths.R")

analysis_dir <- file.path(
  picrust2_results_dir,
  "maaslin3",
  "pathways_primary"
)

results_file <- file.path(
  analysis_dir,
  "pathway_group_results_with_descriptions.tsv"
)

output_dir <- file.path(
  analysis_dir,
  "custom_figures"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Number of unique pathways displayed.
top_n_pathways <- 25

comparison_levels <- c(
  "ABH vs H",
  "AAH vs H"
)

# q >= 0.05 is white.
q_threshold <- 0.05
q_floor <- 0.000001

q_limits <- -log10(
  c(q_threshold, q_floor)
)

q_breaks <- -log10(
  c(
    0.05,
    0.01,
    0.001,
    0.0001,
    0.000001
  )
)

q_labels <- c(
  "≥0.05",
  "0.01",
  "0.001",
  "0.0001",
  "≤0.000001"
)

q_colour_positions <- scales::rescale(
  c(
    q_limits[1],
    2,
    4,
    q_limits[2]
  )
)

abundance_colours <- c(
  "#FFFFFF",
  "#E7CEE8",
  "#C06FBE",
  "#8E138A"
)

prevalence_colours <- c(
  "#FFFFFF",
  "#D7ECEA",
  "#69B8B2",
  "#168C87"
)

q_colour_value <- function(q) {

  ifelse(
    !is.na(q) & q < q_threshold,
    -log10(pmax(q, q_floor)),
    q_limits[1]
  )
}

# Use the correct minus sign in axis labels.
minus_labels <- function(x) {

  labels <- scales::label_number(
    accuracy = 0.1,
    trim = TRUE
  )(x)

  sub("^-", "−", labels)
}

# Replace HTML entities in MetaCyc descriptions.
clean_pathway_description <- function(x) {

  x |>
    str_replace_all(
      c(
        "&alpha;" = "α",
        "&beta;" = "β",
        "&gamma;" = "γ",
        "&delta;" = "δ",
        "&amp;" = "&"
      )
    ) |>
    str_squish()
}

# Read MaAsLin3 results.
results <- read_tsv(
  results_file,
  show_col_types = FALSE
) |>
  mutate(
    description = clean_pathway_description(
      description
    ),
    comparison = factor(
      comparison,
      levels = comparison_levels
    ),
    model = factor(
      model,
      levels = c(
        "abundance",
        "prevalence"
      )
    )
  )

if (!all(
  comparison_levels %in%
    as.character(results$comparison)
)) {
  stop(
    "Expected group comparisons were not found."
  )
}

# Select the 25 pathways with the smallest
# individual q-values.
selected_pathways <- results |>
  filter(
    !is.na(qval_individual),
    qval_individual < q_threshold
  ) |>
  group_by(
    feature,
    description
  ) |>
  summarise(
    minimum_q = min(
      qval_individual,
      na.rm = TRUE
    ),
    maximum_effect_from_null = max(
      abs(coef - null_hypothesis),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  arrange(
    minimum_q,
    desc(maximum_effect_from_null),
    feature
  ) |>
  slice_head(
    n = top_n_pathways
  ) |>
  mutate(
    selection_rank = row_number(),
    pathway_label = paste0(
      description,
      " (",
      feature,
      ")"
    ) |>
      str_wrap(width = 52)
  )

if (nrow(selected_pathways) == 0) {
  stop(
    "No significant group-associated pathways found."
  )
}

# Save the pathway-selection table.
write_tsv(
  selected_pathways |>
    select(
      selection_rank,
      feature,
      description,
      minimum_q,
      maximum_effect_from_null
    ),
  file.path(
    output_dir,
    "top25_pathway_selection.tsv"
  )
)

pathway_order <- selected_pathways$pathway_label

# Prepare plotting data.
plot_data <- results |>
  filter(
    feature %in% selected_pathways$feature
  ) |>
  left_join(
    selected_pathways |>
      select(
        feature,
        pathway_label
      ),
    by = "feature"
  ) |>
  mutate(
    pathway_label = factor(
      pathway_label,
      levels = rev(pathway_order)
    ),
    confidence_low = coef - 1.96 * stderr,
    confidence_high = coef + 1.96 * stderr,
    q_colour = q_colour_value(
      qval_individual
    )
  )

valid_data <- plot_data |>
  filter(
    is.finite(coef),
    is.finite(stderr),
    is.finite(confidence_low),
    is.finite(confidence_high)
  )

# Abundance is tested against a median-based null.
# Prevalence is tested against zero.
null_lines <- plot_data |>
  filter(
    is.finite(null_hypothesis)
  ) |>
  distinct(
    comparison,
    model,
    null_hypothesis
  )

# Use symmetrical x-axis limits.
x_values <- c(
  valid_data$confidence_low,
  valid_data$confidence_high,
  null_lines$null_hypothesis
)

x_limit <- ceiling(
  max(
    abs(x_values),
    na.rm = TRUE
  ) * 1.08
)

if (!is.finite(x_limit) || x_limit < 1) {
  x_limit <- 1
}

# Save exactly the data displayed in the figure.
write_tsv(
  plot_data |>
    mutate(
      pathway_label = as.character(
        pathway_label
      ),
      comparison = as.character(
        comparison
      ),
      model = as.character(
        model
      )
    ),
  file.path(
    output_dir,
    "top25_pathway_plot_data.tsv"
  )
)

# Create figure.
pathway_plot <- ggplot() +

  # Null-hypothesis lines.
  geom_vline(
    data = null_lines,
    aes(
      xintercept = null_hypothesis,
      linetype = model,
      colour = model
    ),
    linewidth = 0.65
  ) +

  # 95% confidence intervals.
  geom_errorbar(
    data = valid_data,
    aes(
      y = pathway_label,
      xmin = confidence_low,
      xmax = confidence_high
    ),
    orientation = "y",
    width = 0.16,
    linewidth = 0.55,
    colour = "black"
  ) +

  # Abundance associations.
  geom_point(
    data = filter(
      valid_data,
      model == "abundance"
    ),
    aes(
      x = coef,
      y = pathway_label,
      fill = q_colour,
      shape = model
    ),
    colour = "black",
    stroke = 0.7,
    size = 3.4
  ) +

  scale_fill_gradientn(
    name = "Abundance q",
    colours = abundance_colours,
    values = q_colour_positions,
    limits = q_limits,
    breaks = q_breaks,
    labels = q_labels,
    oob = scales::squish,
    na.value = "white",
    guide = guide_colorbar(
      order = 3,
      barheight = grid::unit(
        42,
        "mm"
      ),
      barwidth = grid::unit(
        5.5,
        "mm"
      ),
      title.position = "top"
    )
  ) +

  ggnewscale::new_scale_fill() +

  # Prevalence associations.
  geom_point(
    data = filter(
      valid_data,
      model == "prevalence"
    ),
    aes(
      x = coef,
      y = pathway_label,
      fill = q_colour,
      shape = model
    ),
    colour = "black",
    stroke = 0.7,
    size = 3.6
  ) +

  scale_fill_gradientn(
    name = "Prevalence q",
    colours = prevalence_colours,
    values = q_colour_positions,
    limits = q_limits,
    breaks = q_breaks,
    labels = q_labels,
    oob = scales::squish,
    na.value = "white",
    guide = guide_colorbar(
      order = 4,
      barheight = grid::unit(
        42,
        "mm"
      ),
      barwidth = grid::unit(
        5.5,
        "mm"
      ),
      title.position = "top"
    )
  ) +

  scale_shape_manual(
    name = "Association",
    values = c(
      abundance = 21,
      prevalence = 24
    ),
    labels = c(
      abundance = "Abundance",
      prevalence = "Prevalence"
    ),
    guide = guide_legend(
      order = 2,
      override.aes = list(
        fill = "white",
        size = 3.5
      )
    )
  ) +

  scale_linetype_manual(
    name = "Null hypothesis",
    values = c(
      abundance = "dashed",
      prevalence = "solid"
    ),
    labels = c(
      abundance = "Abundance",
      prevalence = "Prevalence"
    ),
    guide = guide_legend(
      order = 1
    )
  ) +

  scale_colour_manual(
    name = "Null hypothesis",
    values = c(
      abundance = "grey55",
      prevalence = "black"
    ),
    labels = c(
      abundance = "Abundance",
      prevalence = "Prevalence"
    ),
    guide = guide_legend(
      order = 1
    )
  ) +

  scale_x_continuous(
    limits = c(
      -x_limit,
      x_limit
    ),
    labels = minus_labels,
    expand = expansion(
      mult = 0.03
    )
  ) +

  facet_grid(
    . ~ comparison,
    drop = FALSE
  ) +

  labs(
    x = "β coefficient",
    y = "Predicted MetaCyc pathway"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    panel.grid.minor = element_blank(),

    panel.grid.major.y = element_line(
      colour = "grey90",
      linewidth = 0.35
    ),

    panel.spacing.x = grid::unit(
      3,
      "mm"
    ),

    strip.background = element_rect(
      fill = "white",
      colour = "black"
    ),

    strip.text = element_text(
      size = 15,
      face = "bold"
    ),

    axis.text.x = element_text(
      size = 12,
      colour = "black"
    ),

    axis.text.y = element_text(
      size = 10.5,
      colour = "black",
      lineheight = 0.95
    ),

    axis.title = element_text(
      size = 15
    ),

    legend.position = "right",

    legend.title = element_text(
      size = 14,
      face = "bold"
    ),

    legend.text = element_text(
      size = 12
    ),

    legend.box = "vertical",

    plot.margin = margin(
      10,
      10,
      10,
      10
    )
  )

# Save PNG.
ggsave(
  file.path(
    output_dir,
    "maaslin3_predicted_pathways_top25.png"
  ),
  pathway_plot,
  width = 18,
  height = 14.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

# Save vector PDF.
ggsave(
  file.path(
    output_dir,
    "maaslin3_predicted_pathways_top25.pdf"
  ),
  pathway_plot,
  device = cairo_pdf,
  width = 18,
  height = 14.5,
  units = "in",
  bg = "white"
)

# Save a short visualization summary.
write_tsv(
  tibble::tibble(
    parameter = c(
      "Unique significant pathways before selection",
      "Pathways displayed",
      "Selection criterion",
      "Individual q threshold",
      "Smallest displayed q"
    ),
    value = c(
      results |>
        filter(
          qval_individual < q_threshold
        ) |>
        distinct(feature) |>
        nrow(),

      nrow(selected_pathways),

      paste(
        "Minimum individual q across",
        "both comparisons and models"
      ),

      q_threshold,

      min(selected_pathways$minimum_q)
    )
  ),
  file.path(
    output_dir,
    "custom_figure_summary.tsv"
  )
)

message(
  "Custom pathway MaAsLin3 figure completed."
)

message(
  "Pathways displayed: ",
  nrow(selected_pathways)
)

message(
  "Output directory: ",
  output_dir
)
