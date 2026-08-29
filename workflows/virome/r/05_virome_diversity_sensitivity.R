# Sensitivity analyses of family-level virome diversity.
#
# Four datasets are compared:
# 1. All samples, including Parvoviridae.
# 2. Without the low-depth sample RNA_S12046Nr12.
# 3. Without Parvoviridae.
# 4. Without both RNA_S12046Nr12 and Parvoviridae.
#
# Rarefaction is not performed.

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
  "sensitivity"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

minimum_total_reads <- 10
low_depth_sample <- "RNA_S12046Nr12"
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
)

family_summary <- read_tsv(
  file.path(
    prepared_dir,
    "viral_family_summary.tsv"
  ),
  show_col_types = FALSE
)

base_families <- family_summary |>
  filter(
    family != "Unclassified viruses",
    total_reads >= minimum_total_reads
  ) |>
  pull(family)

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

run_sensitivity_analysis <- function(
  analysis_id,
  analysis_label,
  remove_nr12,
  remove_parvoviridae
) {

  analysis_metadata <- metadata

  if (remove_nr12) {
    analysis_metadata <- analysis_metadata |>
      filter(
        sample_id != low_depth_sample
      )
  }

  retained_families <- base_families

  if (remove_parvoviridae) {
    retained_families <- setdiff(
      retained_families,
      "Parvoviridae"
    )
  }

  analysis_counts <- counts |>
    filter(
      sample_id %in%
        analysis_metadata$sample_id
    ) |>
    arrange(
      match(
        sample_id,
        analysis_metadata$sample_id
      )
    )

  count_matrix <- analysis_counts |>
    select(
      all_of(retained_families)
    ) |>
    as.matrix()

  rownames(count_matrix) <-
    analysis_counts$sample_id

  storage.mode(count_matrix) <- "numeric"

  relative_matrix <-
    count_matrix / rowSums(count_matrix)

  hellinger_matrix <- sqrt(relative_matrix)

  # Alpha-diversity.
  shannon <- vegan::diversity(
    relative_matrix,
    index = "shannon"
  )

  simpson <- vegan::diversity(
    relative_matrix,
    index = "simpson"
  )

  observed_families <- vegan::specnumber(
    relative_matrix
  )

  pielou <- ifelse(
    observed_families > 1,
    shannon / log(observed_families),
    NA_real_
  )

  alpha_values <- tibble(
    analysis_id = analysis_id,
    analysis_label = analysis_label,
    sample_id = rownames(relative_matrix),
    shannon = as.numeric(shannon),
    simpson = as.numeric(simpson),
    pielou_evenness = as.numeric(pielou),
    observed_families =
      as.numeric(observed_families)
  ) |>
    left_join(
      analysis_metadata,
      by = "sample_id"
    )

  alpha_long <- alpha_values |>
    pivot_longer(
      cols = c(
        shannon,
        simpson,
        pielou_evenness
      ),
      names_to = "metric",
      values_to = "value"
    )

  alpha_kruskal <- bind_rows(
    lapply(
      c(
        "shannon",
        "simpson",
        "pielou_evenness"
      ),
      function(metric_name) {

        metric_data <- alpha_long |>
          filter(
            metric == metric_name
          )

        test <- kruskal.test(
          value ~ group,
          data = metric_data
        )

        tibble(
          analysis_id = analysis_id,
          analysis_label = analysis_label,
          metric = metric_name,
          statistic =
            unname(test$statistic),
          degrees_of_freedom =
            unname(test$parameter),
          p_value = test$p.value
        )
      }
    )
  )

  alpha_depth <- bind_rows(
    lapply(
      c(
        "shannon",
        "simpson",
        "pielou_evenness"
      ),
      function(metric_name) {

        metric_values <-
          alpha_values[[metric_name]]

        test <- cor.test(
          metric_values,
          log10(
            alpha_values$viral_pairs + 1
          ),
          method = "spearman",
          exact = FALSE
        )

        tibble(
          analysis_id = analysis_id,
          analysis_label = analysis_label,
          metric = metric_name,
          rho = unname(test$estimate),
          p_value = test$p.value
        )
      }
    )
  )

  # Beta-diversity.
  bray_curtis <- vegan::vegdist(
    hellinger_matrix,
    method = "bray"
  )

  pcoa <- vegan::wcmdscale(
    bray_curtis,
    k = 2,
    eig = TRUE,
    add = "lingoes"
  )

  positive_eigenvalues <- pcoa$eig[
    pcoa$eig > 0
  ]

  axis_percent <- 100 *
    pcoa$eig[1:2] /
    sum(positive_eigenvalues)

  pcoa_data <- tibble(
    analysis_id = analysis_id,
    analysis_label = analysis_label,
    sample_id = rownames(pcoa$points),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2],
    axis_1_percent = axis_percent[1],
    axis_2_percent = axis_percent[2]
  ) |>
    left_join(
      analysis_metadata,
      by = "sample_id"
    )

  set.seed(random_seed)

  permanova_fit <- vegan::adonis2(
    bray_curtis ~ group,
    data = analysis_metadata,
    permutations = permutations
  )

  permanova <- tibble(
    analysis_id = analysis_id,
    analysis_label = analysis_label,
    samples = nrow(analysis_metadata),
    families = length(retained_families),
    r_squared = permanova_fit$R2[1],
    pseudo_f = permanova_fit$F[1],
    p_value =
      permanova_fit$`Pr(>F)`[1]
  )

  dispersion_fit <- vegan::betadisper(
    bray_curtis,
    analysis_metadata$group,
    type = "median",
    add = "lingoes"
  )

  set.seed(random_seed)

  dispersion_test <- permutest(
    dispersion_fit,
    permutations = permutations
  )

  permdisp <- tibble(
    analysis_id = analysis_id,
    analysis_label = analysis_label,
    f_statistic =
      dispersion_test$tab$F[1],
    p_value =
      dispersion_test$tab$`Pr(>F)`[1]
  )

  settings <- tibble(
    analysis_id = analysis_id,
    analysis_label = analysis_label,
    samples = nrow(analysis_metadata),
    hibernation_samples = sum(
      analysis_metadata$group ==
        "hibernation"
    ),
    families = length(retained_families),
    nr12_included = !remove_nr12,
    parvoviridae_included =
      !remove_parvoviridae
  )

  list(
    settings = settings,
    alpha_values = alpha_values,
    alpha_kruskal = alpha_kruskal,
    alpha_depth = alpha_depth,
    pcoa = pcoa_data,
    permanova = permanova,
    permdisp = permdisp
  )
}

analysis_settings <- list(
  list(
    id = "primary",
    label = "All samples",
    remove_nr12 = FALSE,
    remove_parvoviridae = FALSE
  ),
  list(
    id = "without_nr12",
    label = "Without H2 (Nr12)",
    remove_nr12 = TRUE,
    remove_parvoviridae = FALSE
  ),
  list(
    id = "without_parvoviridae",
    label = "Without Parvoviridae",
    remove_nr12 = FALSE,
    remove_parvoviridae = TRUE
  ),
  list(
    id = "without_nr12_and_parvoviridae",
    label = "Without H2 and Parvoviridae",
    remove_nr12 = TRUE,
    remove_parvoviridae = TRUE
  )
)

results <- lapply(
  analysis_settings,
  function(setting) {

    run_sensitivity_analysis(
      analysis_id = setting$id,
      analysis_label = setting$label,
      remove_nr12 = setting$remove_nr12,
      remove_parvoviridae =
        setting$remove_parvoviridae
    )
  }
)

settings_table <- bind_rows(
  lapply(results, `[[`, "settings")
)

alpha_values_table <- bind_rows(
  lapply(results, `[[`, "alpha_values")
)

alpha_kruskal_table <- bind_rows(
  lapply(results, `[[`, "alpha_kruskal")
)

alpha_depth_table <- bind_rows(
  lapply(results, `[[`, "alpha_depth")
)

pcoa_table <- bind_rows(
  lapply(results, `[[`, "pcoa")
)

permanova_table <- bind_rows(
  lapply(results, `[[`, "permanova")
)

permdisp_table <- bind_rows(
  lapply(results, `[[`, "permdisp")
)

write_tsv(
  settings_table,
  file.path(
    output_dir,
    "sensitivity_settings.tsv"
  )
)

write_tsv(
  alpha_values_table |>
    mutate(group = as.character(group)),
  file.path(
    output_dir,
    "sensitivity_alpha_values.tsv"
  )
)

write_tsv(
  alpha_kruskal_table,
  file.path(
    output_dir,
    "sensitivity_alpha_kruskal_wallis.tsv"
  )
)

write_tsv(
  alpha_depth_table,
  file.path(
    output_dir,
    "sensitivity_alpha_depth_correlations.tsv"
  )
)

write_tsv(
  pcoa_table |>
    mutate(group = as.character(group)),
  file.path(
    output_dir,
    "sensitivity_pcoa_coordinates.tsv"
  )
)

write_tsv(
  permanova_table,
  file.path(
    output_dir,
    "sensitivity_permanova.tsv"
  )
)

write_tsv(
  permdisp_table,
  file.path(
    output_dir,
    "sensitivity_permdisp.tsv"
  )
)

# Create one combined statistical summary.
alpha_summary_wide <- alpha_kruskal_table |>
  select(
    analysis_id,
    metric,
    p_value
  ) |>
  pivot_wider(
    names_from = metric,
    values_from = p_value,
    names_prefix = "kw_p_"
  )

combined_summary <- settings_table |>
  left_join(
    alpha_summary_wide,
    by = "analysis_id"
  ) |>
  left_join(
    permanova_table |>
      select(
        analysis_id,
        permanova_r_squared = r_squared,
        permanova_p_value = p_value
      ),
    by = "analysis_id"
  ) |>
  left_join(
    permdisp_table |>
      select(
        analysis_id,
        permdisp_p_value = p_value
      ),
    by = "analysis_id"
  )

write_tsv(
  combined_summary,
  file.path(
    output_dir,
    "sensitivity_diversity_summary.tsv"
  )
)

# ============================================================
# Diagnostic PCoA figure
# ============================================================

make_pcoa_plot <- function(result) {

  plot_data <- result$pcoa
  permanova_result <- result$permanova
  permdisp_result <- result$permdisp

  axis_1 <- unique(
    plot_data$axis_1_percent
  )

  axis_2 <- unique(
    plot_data$axis_2_percent
  )

  statistic_label <- paste0(
    "PERMANOVA: R² = ",
    formatC(
      permanova_result$r_squared,
      digits = 3,
      format = "f"
    ),
    ", ",
    format_p(
      permanova_result$p_value
    ),
    "\nPERMDISP: ",
    format_p(
      permdisp_result$p_value
    )
  )

  ggplot(
    plot_data,
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
      label = statistic_label,
      hjust = -0.03,
      vjust = 1.05,
      size = 3.7,
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
        mult = c(0.05, 0.35)
      )
    ) +
    coord_fixed(
      ratio = 1
    ) +
    labs(
      title = unique(
        plot_data$analysis_label
      ),
      x = paste0(
        "PCoA axis 1 (",
        round(axis_1, 1),
        "%)"
      ),
      y = paste0(
        "PCoA axis 2 (",
        round(axis_2, 1),
        "%)"
      )
    ) +
    theme_classic(
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}

pcoa_plots <- lapply(
  results,
  make_pcoa_plot
)

sensitivity_figure <- wrap_plots(
  pcoa_plots,
  ncol = 2,
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
    "virome_diversity_sensitivity_pcoa.png"
  ),
  sensitivity_figure,
  width = 13,
  height = 10,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "virome_diversity_sensitivity_pcoa.pdf"
  ),
  sensitivity_figure,
  width = 13,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

cat("\nVirome diversity sensitivity analysis completed.\n")
cat("Rarefaction: not performed\n")
cat("Results:", output_dir, "\n\n")

print(combined_summary)
cat("\n")
print(alpha_depth_table)
