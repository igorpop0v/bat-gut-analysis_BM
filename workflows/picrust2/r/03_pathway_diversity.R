# Alpha- and beta-diversity of predicted MetaCyc pathways.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(vegan)

source("workflows/picrust2/r/local_paths.R")

pathway_file <- file.path(
  picrust2_input_dir,
  "metacyc_pathway_abundance_with_descriptions.tsv.gz"
)

metadata_file <- file.path(
  picrust2_input_dir,
  "sample_metadata.tsv"
)

output_dir <- file.path(
  picrust2_results_dir,
  "pathway_diversity"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

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

permutations <- 9999
random_seed <- 10000

# Read metadata.
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

if (any(is.na(metadata$group))) {
  stop("Unexpected or missing group values in metadata.")
}

# Read predicted MetaCyc pathway abundances.
pathway_table <- read_tsv(
  pathway_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

if (anyDuplicated(pathway_table$pathway)) {
  stop("Duplicated MetaCyc pathway IDs were found.")
}

if (!all(metadata$sample_id %in% names(pathway_table))) {
  stop("Some metadata samples are absent from the pathway table.")
}

# Create a sample-by-pathway abundance matrix
# in the same order as the metadata.
pathway_abundance <- t(
  as.matrix(
    pathway_table[, metadata$sample_id]
  )
)

rownames(pathway_abundance) <- metadata$sample_id
colnames(pathway_abundance) <- pathway_table$pathway

if (any(pathway_abundance < 0)) {
  stop("Negative pathway abundances were found.")
}

sample_totals <- rowSums(pathway_abundance)

if (any(sample_totals <= 0)) {
  stop("At least one sample has zero total pathway abundance.")
}

# TSS normalization: the pathway abundances
# in every sample sum to one.
pathway_relative <- pathway_abundance / sample_totals

if (!all(abs(rowSums(pathway_relative) - 1) < 1e-10)) {
  stop("TSS normalization failed.")
}

# Save the complete normalized pathway table.
relative_output <- as.data.frame(
  pathway_relative,
  check.names = FALSE
) |>
  tibble::rownames_to_column("sample_id")

write_tsv(
  relative_output,
  file.path(
    output_dir,
    "pathway_relative_abundance_tss.tsv"
  )
)

# Calculate alpha-diversity using all predicted pathways.
observed_pathways <- vegan::specnumber(
  pathway_relative
)

shannon <- vegan::diversity(
  pathway_relative,
  index = "shannon"
)

simpson <- vegan::diversity(
  pathway_relative,
  index = "simpson"
)

pielou_evenness <- ifelse(
  observed_pathways > 1,
  shannon / log(observed_pathways),
  NA_real_
)

alpha_wide <- tibble(
  sample_id = rownames(pathway_relative),
  shannon = as.numeric(shannon),
  simpson = as.numeric(simpson),
  pielou_evenness = as.numeric(pielou_evenness)
) |>
  left_join(
    metadata,
    by = "sample_id"
  )

write_tsv(
  alpha_wide |>
    mutate(group = as.character(group)),
  file.path(
    output_dir,
    "pathway_alpha_diversity.tsv"
  )
)

alpha_long <- alpha_wide |>
  pivot_longer(
    cols = c(
      shannon,
      simpson,
      pielou_evenness
    ),
    names_to = "metric",
    values_to = "value"
  )

# Overall Kruskal-Wallis tests.
alpha_overall <- bind_rows(
  lapply(
    unique(alpha_long$metric),
    function(metric_name) {

      metric_data <- alpha_long |>
        filter(metric == metric_name)

      test <- kruskal.test(
        value ~ group,
        data = metric_data
      )

      tibble(
        metric = metric_name,
        statistic = unname(test$statistic),
        degrees_of_freedom = unname(test$parameter),
        p_value = test$p.value
      )
    }
  )
)

write_tsv(
  alpha_overall,
  file.path(
    output_dir,
    "pathway_alpha_kruskal_wallis.tsv"
  )
)

# All three pairwise group comparisons.
group_pairs <- combn(
  group_levels,
  2,
  simplify = FALSE
)

# Pairwise Wilcoxon tests.
# BH correction is applied independently within each metric.
alpha_pairwise <- bind_rows(
  lapply(
    unique(alpha_long$metric),
    function(metric_name) {

      metric_data <- alpha_long |>
        filter(metric == metric_name)

      pair_results <- bind_rows(
        lapply(
          group_pairs,
          function(pair) {

            first_values <- metric_data |>
              filter(group == pair[1]) |>
              pull(value)

            second_values <- metric_data |>
              filter(group == pair[2]) |>
              pull(value)

            test <- wilcox.test(
              first_values,
              second_values,
              exact = FALSE
            )

            tibble(
              metric = metric_name,
              group_1 = pair[1],
              group_2 = pair[2],
              statistic = unname(test$statistic),
              p_value = test$p.value
            )
          }
        )
      )

      pair_results |>
        mutate(
          q_value = p.adjust(
            p_value,
            method = "BH"
          )
        )
    }
  )
)

write_tsv(
  alpha_pairwise,
  file.path(
    output_dir,
    "pathway_alpha_pairwise_wilcoxon.tsv"
  )
)

# Bray-Curtis distance matrix using all pathways.
bray_curtis <- vegan::vegdist(
  pathway_relative,
  method = "bray"
)

bray_matrix <- as.data.frame(
  as.matrix(bray_curtis),
  check.names = FALSE
) |>
  tibble::rownames_to_column("sample_id")

write_tsv(
  bray_matrix,
  file.path(
    output_dir,
    "pathway_bray_curtis_distance_matrix.tsv"
  )
)

# PCoA with a Lingoes correction for negative eigenvalues.
pcoa <- vegan::wcmdscale(
  bray_curtis,
  k = 2,
  eig = TRUE,
  add = "lingoes"
)

positive_eigenvalues <- pcoa$eig[
  pcoa$eig > 0
]

axis_percent <- 100 * pcoa$eig[1:2] /
  sum(positive_eigenvalues)

pcoa_data <- tibble(
  sample_id = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2]
) |>
  left_join(
    metadata,
    by = "sample_id"
  )

write_tsv(
  pcoa_data |>
    mutate(group = as.character(group)),
  file.path(
    output_dir,
    "pathway_bray_curtis_pcoa.tsv"
  )
)

# Overall PERMANOVA.
set.seed(random_seed)

permanova_fit <- vegan::adonis2(
  bray_curtis ~ group,
  data = metadata,
  permutations = permutations
)

permanova_overall <- tibble(
  method = "PERMANOVA",
  permutations = permutations,
  degrees_of_freedom = permanova_fit$Df[1],
  sum_of_squares = permanova_fit$SumOfSqs[1],
  r_squared = permanova_fit$R2[1],
  pseudo_f = permanova_fit$F[1],
  p_value = permanova_fit$`Pr(>F)`[1]
)

write_tsv(
  permanova_overall,
  file.path(
    output_dir,
    "pathway_bray_curtis_permanova.tsv"
  )
)

# Pairwise PERMANOVA.
permanova_pairwise <- bind_rows(
  lapply(
    seq_along(group_pairs),
    function(pair_number) {

      pair <- group_pairs[[pair_number]]

      pair_metadata <- metadata |>
        filter(group %in% pair) |>
        mutate(
          group = droplevels(group)
        )

      pair_distance <- as.dist(
        as.matrix(bray_curtis)[
          pair_metadata$sample_id,
          pair_metadata$sample_id
        ]
      )

      set.seed(
        random_seed + pair_number
      )

      pair_fit <- vegan::adonis2(
        pair_distance ~ group,
        data = pair_metadata,
        permutations = permutations
      )

      tibble(
        group_1 = pair[1],
        group_2 = pair[2],
        degrees_of_freedom = pair_fit$Df[1],
        r_squared = pair_fit$R2[1],
        pseudo_f = pair_fit$F[1],
        p_value = pair_fit$`Pr(>F)`[1]
      )
    }
  )
) |>
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    )
  )

write_tsv(
  permanova_pairwise,
  file.path(
    output_dir,
    "pathway_bray_curtis_pairwise_permanova.tsv"
  )
)

# Overall PERMDISP using distances to group medians.
dispersion_fit <- vegan::betadisper(
  bray_curtis,
  metadata$group,
  type = "median",
  add = "lingoes"
)

set.seed(random_seed)

dispersion_test <- permutest(
  dispersion_fit,
  permutations = permutations
)

permdisp_overall <- tibble(
  method = "PERMDISP",
  permutations = permutations,
  degrees_of_freedom = dispersion_test$tab$Df[1],
  f_statistic = dispersion_test$tab$F[1],
  p_value = dispersion_test$tab$`Pr(>F)`[1]
)

write_tsv(
  permdisp_overall,
  file.path(
    output_dir,
    "pathway_bray_curtis_permdisp.tsv"
  )
)

# Pairwise PERMDISP.
permdisp_pairwise <- bind_rows(
  lapply(
    seq_along(group_pairs),
    function(pair_number) {

      pair <- group_pairs[[pair_number]]

      pair_metadata <- metadata |>
        filter(group %in% pair) |>
        mutate(
          group = droplevels(group)
        )

      pair_distance <- as.dist(
        as.matrix(bray_curtis)[
          pair_metadata$sample_id,
          pair_metadata$sample_id
        ]
      )

      pair_dispersion <- vegan::betadisper(
        pair_distance,
        pair_metadata$group,
        type = "median",
        add = "lingoes"
      )

      set.seed(
        random_seed + 100 + pair_number
      )

      pair_test <- permutest(
        pair_dispersion,
        permutations = permutations
      )

      tibble(
        group_1 = pair[1],
        group_2 = pair[2],
        degrees_of_freedom = pair_test$tab$Df[1],
        f_statistic = pair_test$tab$F[1],
        p_value = pair_test$tab$`Pr(>F)`[1]
      )
    }
  )
) |>
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    )
  )

write_tsv(
  permdisp_pairwise,
  file.path(
    output_dir,
    "pathway_bray_curtis_pairwise_permdisp.tsv"
  )
)

# Formatting functions.
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

minus_labels <- function(x) {

  sub(
    "^-",
    "−",
    scales::label_number(
      trim = TRUE
    )(x)
  )
}

# Create one alpha-diversity panel.
make_alpha_plot <- function(metric_name, y_title) {

  plot_data <- alpha_long |>
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
    ) |>
    arrange(
      x1,
      x2
    )

  y_min <- min(
    plot_data$value,
    na.rm = TRUE
  )

  y_max <- max(
    plot_data$value,
    na.rm = TRUE
  )

  y_range <- y_max - y_min

  if (y_range == 0) {
    y_range <- max(
      abs(y_max),
      1
    ) * 0.05
  }

  # Put significant pairwise comparisons
  # at different vertical positions.
  comparisons <- comparisons |>
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
      x = group,
      y = value,
      fill = group
    )
  ) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.55,
      linewidth = 0.8,
      color = "black",
      show.legend = FALSE
    ) +
    geom_point(
      shape = 21,
      color = "black",
      size = 2.5,
      stroke = 0.55,
      position = position_jitter(
        width = 0.22,
        height = 0,
        seed = random_seed
      ),
      show.legend = FALSE
    ) +
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
    annotate(
      "text",
      x = 1,
      y = overall_y,
      label = paste0(
        "K-W, ",
        format_p(
          overall$p_value
        )
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

# Prepare a clear PCoA statistical annotation.

pairwise_display <- permanova_pairwise |>
  mutate(
    comparison = case_when(
      group_1 == "active" &
        group_2 == "hibernation" ~
        "ABH vs. H",
      
      group_1 == "hibernation" &
        group_2 == "post_hibernation" ~
        "H vs. AAH",
      
      group_1 == "active" &
        group_2 == "post_hibernation" ~
        "ABH vs. AAH",
      
      TRUE ~ paste(
        group_1,
        "vs.",
        group_2
      )
    ),
    
    comparison_order = case_when(
      comparison == "ABH vs. H" ~ 1,
      comparison == "H vs. AAH" ~ 2,
      comparison == "ABH vs. AAH" ~ 3,
      TRUE ~ 4
    ),
    
    q_text = vapply(
      q_value,
      format_q,
      character(1)
    )
  ) |>
  arrange(comparison_order) |>
  mutate(
    result_line = paste0(
      comparison,
      ", ",
      q_text
    )
  ) |>
  pull(result_line) |>
  paste(
    collapse = "\n"
  )

beta_label <- paste0(
  "PERMANOVA: pseudo-F = ",
  formatC(
    permanova_overall$pseudo_f,
    digits = 2,
    format = "f"
  ),
  ", R² = ",
  formatC(
    permanova_overall$r_squared,
    digits = 3,
    format = "f"
  ),
  ", ",
  format_p(
    permanova_overall$p_value
  ),
  "\nPairwise PERMANOVA:\n",
  pairwise_display,
  "\nPERMDISP: F = ",
  formatC(
    permdisp_overall$f_statistic,
    digits = 2,
    format = "f"
  ),
  ", ",
  format_p(
    permdisp_overall$p_value
  )
)

# Bray-Curtis PCoA.
plot_pcoa <- ggplot(
  pcoa_data,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = group
  )
) +
  geom_point(
    size = 2.5
  ) +
  stat_ellipse(
    geom = "path",
    level = 0.95,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
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
  scale_x_continuous(
    labels = minus_labels,
    expand = expansion(
      mult = c(
        0.05,
        0.05
      )
    )
  ) +
  scale_y_continuous(
    labels = minus_labels,
    expand = expansion(
      mult = c(
        0.05,
        0.55
      )
    )
  ) +
  labs(
    x = paste0(
      "PCoA axis 1 (",
      round(
        axis_percent[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PCoA axis 2 (",
      round(
        axis_percent[2],
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
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 16,
        size = 3,
        linewidth = 0
      )
    )
  )

# Combine alpha and beta panels.
final_figure <- (
  plot_shannon |
    plot_simpson |
    plot_pielou
) /
  plot_pcoa +
  plot_layout(
    heights = c(
      1,
      1.15
    ),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(
    output_dir,
    "predicted_pathway_alpha_beta_diversity.png"
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
    "predicted_pathway_alpha_beta_diversity.pdf"
  ),
  final_figure,
  width = 15,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

message("")
message("Predicted pathway diversity analysis completed.")
message("Results: ", output_dir)
message("")
print(alpha_overall)
message("")
print(alpha_pairwise)
message("")
print(permanova_overall)
message("")
print(permanova_pairwise)
message("")
print(permdisp_overall)
message("")
print(permdisp_pairwise)
