# Primary family-level virome diversity analysis.
#
# All 10 samples are retained.
# Rarefaction is not performed.
# Unclassified viruses and extremely low-support families
# are excluded before diversity calculation.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(vegan)

source("workflows/virome/r/local_paths.R")

prepared_dir <- file.path(
  virome_results_dir,
  "prepared_tables"
)

output_dir <- file.path(
  virome_results_dir,
  "diversity",
  "primary_all_samples"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

minimum_total_reads <- 10
permutations <- 9999
random_seed <- 10000

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

# ============================================================
# Read inputs
# ============================================================

metadata <- read_tsv(
  file.path(
    prepared_dir,
    "sample_metadata_and_qc.tsv"
  ),
  show_col_types = FALSE
) |>
  mutate(
    group = factor(
      group,
      levels = group_levels
    )
  )

counts <- read_tsv(
  file.path(
    prepared_dir,
    "viral_family_counts_all_samples.tsv"
  ),
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  arrange(
    match(
      sample_id,
      metadata$sample_id
    )
  )

family_summary <- read_tsv(
  file.path(
    prepared_dir,
    "viral_family_summary.tsv"
  ),
  show_col_types = FALSE
)

# Retain classified families with at least 10 reads
# across the complete dataset.
retained_families <- family_summary |>
  filter(
    family != "Unclassified viruses",
    total_reads >= minimum_total_reads
  ) |>
  pull(family)

count_matrix <- counts |>
  select(
    all_of(retained_families)
  ) |>
  as.matrix()

rownames(count_matrix) <- counts$sample_id

storage.mode(count_matrix) <- "numeric"

# TSS normalization.
relative_matrix <- count_matrix / rowSums(count_matrix)

# Hellinger transformation.
hellinger_matrix <- sqrt(relative_matrix)

write_tsv(
  as.data.frame(
    relative_matrix,
    check.names = FALSE
  ) |>
    tibble::rownames_to_column("sample_id"),
  file.path(
    output_dir,
    "viral_family_relative_abundance.tsv"
  )
)

filter_summary <- tibble(
  parameter = c(
    "Samples",
    "Families before filtering",
    "Families after filtering",
    "Minimum total reads",
    "Unclassified viruses included",
    "Rarefaction",
    "Beta-diversity transformation",
    "Beta-diversity distance"
  ),
  value = c(
    nrow(relative_matrix),
    nrow(family_summary),
    length(retained_families),
    minimum_total_reads,
    "No",
    "No",
    "Hellinger",
    "Bray-Curtis"
  )
)

write_tsv(
  filter_summary,
  file.path(
    output_dir,
    "diversity_analysis_summary.tsv"
  )
)

# ============================================================
# Alpha-diversity
# ============================================================

shannon_values <- vegan::diversity(
  relative_matrix,
  index = "shannon"
)

simpson_values <- vegan::diversity(
  relative_matrix,
  index = "simpson"
)

observed_families <- vegan::specnumber(
  relative_matrix
)

pielou_values <- ifelse(
  observed_families > 1,
  shannon_values / log(observed_families),
  NA_real_
)

alpha_diversity <- tibble(
  sample_id = rownames(relative_matrix),
  shannon = as.numeric(shannon_values),
  simpson = as.numeric(simpson_values),
  pielou_evenness = as.numeric(pielou_values),
  observed_families = as.numeric(observed_families)
) |>
  left_join(
    metadata,
    by = "sample_id"
  )

write_tsv(
  alpha_diversity |>
    mutate(group = as.character(group)),
  file.path(
    output_dir,
    "viral_alpha_diversity.tsv"
  )
)

alpha_long <- alpha_diversity |>
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
    c("shannon", "simpson", "pielou_evenness"),
    function(metric_name) {

      test_data <- alpha_long |>
        filter(metric == metric_name)

      test <- kruskal.test(
        value ~ group,
        data = test_data
      )

      tibble(
        metric = metric_name,
        statistic = unname(test$statistic),
        degrees_of_freedom = unname(
          test$parameter
        ),
        p_value = test$p.value
      )
    }
  )
)

write_tsv(
  alpha_overall,
  file.path(
    output_dir,
    "viral_alpha_kruskal_wallis.tsv"
  )
)

group_pairs <- combn(
  group_levels,
  2,
  simplify = FALSE
)

# Pairwise Wilcoxon tests.
alpha_pairwise <- bind_rows(
  lapply(
    c("shannon", "simpson", "pielou_evenness"),
    function(metric_name) {

      metric_data <- alpha_long |>
        filter(metric == metric_name)

      results <- bind_rows(
        lapply(
          group_pairs,
          function(pair) {

            test <- wilcox.test(
              metric_data |>
                filter(group == pair[1]) |>
                pull(value),

              metric_data |>
                filter(group == pair[2]) |>
                pull(value),

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

      results |>
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
    "viral_alpha_pairwise_wilcoxon.tsv"
  )
)

# Association between alpha-diversity and viral depth.
alpha_depth_correlations <- bind_rows(
  lapply(
    c("shannon", "simpson", "pielou_evenness"),
    function(metric_name) {

      metric_values <- alpha_diversity[[metric_name]]

      test <- cor.test(
        metric_values,
        log10(alpha_diversity$viral_pairs + 1),
        method = "spearman",
        exact = FALSE
      )

      tibble(
        metric = metric_name,
        method = "Spearman",
        rho = unname(test$estimate),
        p_value = test$p.value
      )
    }
  )
)

write_tsv(
  alpha_depth_correlations,
  file.path(
    output_dir,
    "viral_alpha_depth_correlations.tsv"
  )
)

# ============================================================
# Beta-diversity
# ============================================================

bray_curtis <- vegan::vegdist(
  hellinger_matrix,
  method = "bray"
)

write_tsv(
  as.data.frame(
    as.matrix(bray_curtis),
    check.names = FALSE
  ) |>
    tibble::rownames_to_column("sample_id"),
  file.path(
    output_dir,
    "viral_hellinger_bray_curtis_matrix.tsv"
  )
)

# PCoA with Lingoes correction.
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
    "viral_hellinger_bray_pcoa.tsv"
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
  r_squared = permanova_fit$R2[1],
  pseudo_f = permanova_fit$F[1],
  p_value = permanova_fit$`Pr(>F)`[1]
)

write_tsv(
  permanova_overall,
  file.path(
    output_dir,
    "viral_permanova.tsv"
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

      fit <- vegan::adonis2(
        pair_distance ~ group,
        data = pair_metadata,
        permutations = permutations
      )

      tibble(
        group_1 = pair[1],
        group_2 = pair[2],
        r_squared = fit$R2[1],
        pseudo_f = fit$F[1],
        p_value = fit$`Pr(>F)`[1]
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
    "viral_pairwise_permanova.tsv"
  )
)

# Overall PERMDISP.
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
  degrees_of_freedom =
    dispersion_test$tab$Df[1],
  f_statistic =
    dispersion_test$tab$F[1],
  p_value =
    dispersion_test$tab$`Pr(>F)`[1]
)

write_tsv(
  permdisp_overall,
  file.path(
    output_dir,
    "viral_permdisp.tsv"
  )
)

# ============================================================
# Visualization
# ============================================================

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

make_alpha_plot <- function(
  metric_name,
  y_title
) {

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
      x1 = match(group_1, group_levels),
      x2 = match(group_2, group_levels)
    )

  y_min <- min(plot_data$value)
  y_max <- max(plot_data$value)
  y_range <- y_max - y_min

  if (y_range == 0) {
    y_range <- 0.05
  }

  comparisons <- comparisons |>
    mutate(
      bracket_y = y_max +
        (
          0.08 +
            0.14 * (row_number() - 1)
        ) * y_range,

      text_y = bracket_y +
        0.05 * y_range
    )

  overall_y <- y_max +
    (
      0.18 +
        0.14 * nrow(comparisons)
    ) * y_range

  upper_limit <- overall_y +
    0.10 * y_range

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
      size = 3,
      stroke = 0.6,
      position = position_jitter(
        width = 0.18,
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
      size = 4.3
    ) +
    annotate(
      "text",
      x = 1,
      y = overall_y,
      label = paste0(
        "K-W, ",
        format_p(overall$p_value)
      ),
      hjust = 0,
      size = 4.5
    ) +
    scale_fill_manual(
      values = group_colors,
      guide = "none"
    ) +
    scale_x_discrete(
      labels = group_labels
    ) +
    coord_cartesian(
      ylim = c(
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
    )
}

shannon_plot <- make_alpha_plot(
  "shannon",
  "Shannon index"
)

simpson_plot <- make_alpha_plot(
  "simpson",
  "Simpson index"
)

pielou_plot <- make_alpha_plot(
  "pielou_evenness",
  "Pielou evenness"
)

if (
  all(
    permanova_pairwise$q_value > 0.05
  )
) {
  
  pairwise_text <-
    "All pairwise comparisons: q > 0.05"
  
} else {
  
  pairwise_text <- permanova_pairwise |>
    mutate(
      comparison = case_when(
        group_1 == "active" &
          group_2 == "hibernation" ~
          "ABH vs. H",
        
        group_1 == "active" &
          group_2 == "post_hibernation" ~
          "ABH vs. AAH",
        
        group_1 == "hibernation" &
          group_2 == "post_hibernation" ~
          "H vs. AAH"
      ),
      
      result = paste0(
        comparison,
        ", ",
        vapply(
          q_value,
          format_q,
          character(1)
        )
      )
    ) |>
    pull(result) |>
    paste(collapse = "\n")
}

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
  format_p(permanova_overall$p_value),
  "\n",
  pairwise_text,
  "\nPERMDISP: F = ",
  formatC(
    permdisp_overall$f_statistic,
    digits = 2,
    format = "f"
  ),
  ", ",
  format_p(permdisp_overall$p_value)
)

pcoa_plot <- ggplot(
  pcoa_data,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = group
  )
) +
  geom_point(
    size = 3
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = beta_label,
    hjust = -0.03,
    vjust = 1.05,
    size = 4,
    fill = "white",
    linewidth = 0
  ) +
  scale_color_manual(
    values = group_colors,
    labels = group_labels,
    name = "Group"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.55)
    )
  ) +
  coord_fixed(
    ratio = 1
  ) +
  labs(
    x = paste0(
      "PCoA axis 1 (",
      round(axis_percent[1], 1),
      "%)"
    ),
    y = paste0(
      "PCoA axis 2 (",
      round(axis_percent[2], 1),
      "%)"
    )
  ) +
  theme_classic(
    base_size = 16
  ) +
  theme(
    legend.title = element_text(
      size = 16
    ),
    legend.text = element_text(
      size = 15
    )
  )

final_figure <- (
  shannon_plot |
    simpson_plot |
    pielou_plot |
    pcoa_plot
) +
  plot_layout(
    widths = c(
      1,
      1,
      1,
      1.35
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
    "primary_virome_alpha_beta_diversity.png"
  ),
  final_figure,
  width = 19,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "primary_virome_alpha_beta_diversity.pdf"
  ),
  final_figure,
  width = 19,
  height = 5.5,
  device = cairo_pdf,
  bg = "white"
)

cat("\nPrimary virome diversity analysis completed.\n")
cat("Rarefaction: not performed\n")
cat("Samples:", nrow(relative_matrix), "\n")
cat(
  "Families retained:",
  length(retained_families),
  "\n"
)
cat("Results:", output_dir, "\n\n")

print(alpha_overall)
cat("\n")
print(alpha_pairwise)
cat("\n")
print(alpha_depth_correlations)
cat("\n")
print(permanova_overall)
cat("\n")
print(permanova_pairwise)
cat("\n")
print(permdisp_overall)
