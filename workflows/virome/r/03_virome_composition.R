# Visualize family-level virome composition.
#
# Outputs:
# 1. Whole-virome composition by sample.
# 2. Mean whole-virome composition by group.
# 3. Descriptive bacteriophage composition in samples with
#    at least 100 phage-classified read pairs.
#
# Group means are calculated from sample-level relative
# abundances, so each sample receives equal weight.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)

source("workflows/virome/r/local_paths.R")

prepared_dir <- file.path(
  virome_results_dir,
  "prepared_tables"
)

analysis_dir <- file.path(
  virome_results_dir,
  "analysis_datasets"
)

output_dir <- file.path(
  virome_results_dir,
  "composition"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

group_levels <- c(
  "ABH",
  "H",
  "AAH"
)

sample_levels <- c(
  "RNA_S12046Nr6",
  "RNA_S12046Nr7",
  "RNA_S12046Nr8",
  "RNA_S12046Nr9",
  "RNA_S12046Nr11",
  "RNA_S12046Nr12",
  "RNA_S12046Nr14",
  "RNA_S12046Nr20",
  "RNA_S12046Nr21",
  "RNA_S12046Nr22"
)

sample_labels <- c(
  "RNA_S12046Nr6" = "ABH1",
  "RNA_S12046Nr7" = "ABH2",
  "RNA_S12046Nr8" = "ABH3",
  "RNA_S12046Nr9" = "ABH4",
  "RNA_S12046Nr11" = "H1",
  "RNA_S12046Nr12" = "H2",
  "RNA_S12046Nr14" = "H3",
  "RNA_S12046Nr20" = "AAH1",
  "RNA_S12046Nr21" = "AAH2",
  "RNA_S12046Nr22" = "AAH3"
)

# Minimum number of phage reads required for a descriptive
# phage-composition bar.
minimum_phage_reads <- 100

metadata <- read_tsv(
  file.path(
    prepared_dir,
    "sample_metadata_and_qc.tsv"
  ),
  show_col_types = FALSE
) |>
  mutate(
    group_label = factor(
      group_label,
      levels = group_levels
    ),
    sample_id = factor(
      sample_id,
      levels = sample_levels
    ),
    sample_label = unname(
      sample_labels[
        as.character(sample_id)
      ]
    )
  )

whole_virome <- read_tsv(
  file.path(
    prepared_dir,
    "viral_family_counts_long_complete.tsv"
  ),
  show_col_types = FALSE
) |>
  left_join(
    metadata |>
      select(
        sample_id,
        sample_label,
        group_label,
        viral_pairs,
        low_depth
      ),
    by = "sample_id"
  )

family_summary <- read_tsv(
  file.path(
    prepared_dir,
    "viral_family_summary.tsv"
  ),
  show_col_types = FALSE
)

# Retain families that are abundant on average or reach at
# least 2% in an individual sample.
displayed_families <- family_summary |>
  filter(
    mean_relative_abundance_pct >= 1 |
      maximum_relative_abundance_pct >= 2
  ) |>
  pull(family)

whole_virome_plot_data <- whole_virome |>
  mutate(
    displayed_family = if_else(
      family %in% displayed_families,
      family,
      "Other"
    )
  ) |>
  group_by(
    sample_id,
    sample_label,
    group_label,
    viral_pairs,
    low_depth,
    displayed_family
  ) |>
  summarise(
    relative_abundance_pct = sum(
      relative_abundance_pct
    ),
    .groups = "drop"
  )

# Order legend by mean relative abundance.
family_order <- whole_virome_plot_data |>
  group_by(displayed_family) |>
  summarise(
    mean_abundance = mean(relative_abundance_pct),
    .groups = "drop"
  ) |>
  arrange(desc(mean_abundance)) |>
  pull(displayed_family)

family_order <- c(
  "Other",
  setdiff(
    family_order,
    "Other"
  )
)

whole_virome_plot_data <- whole_virome_plot_data |>
  mutate(
    displayed_family = factor(
      displayed_family,
      levels = family_order
    )
  )

write_tsv(
  whole_virome_plot_data |>
    mutate(
      sample_id = as.character(sample_id),
      group_label = as.character(group_label),
      displayed_family = as.character(
        displayed_family
      )
    ),
  file.path(
    output_dir,
    "whole_virome_sample_composition.tsv"
  )
)

# Colours for the principal viral families.
whole_virome_colors <- c(
  "Parvoviridae" = "#4D4DFF",
  "Iflaviridae" = "#59A9D3",
  "Autotranscriptaviridae" = "#4E7189",
  "Papillomaviridae" = "#7BA05B",
  "Iridoviridae" = "#F1E487",
  "Pneumoviridae" = "#C16B38",
  "Astroviridae" = "#D33C32",
  "Unclassified viruses" = "#D9D9D9",
  "Other" = "#908A99"
)

# Add a colour for any additional displayed family.
missing_colors <- setdiff(
  family_order,
  names(whole_virome_colors)
)

if (length(missing_colors) > 0) {

  additional_colors <- hue_pal()(
    length(missing_colors)
  )

  names(additional_colors) <- missing_colors

  whole_virome_colors <- c(
    whole_virome_colors,
    additional_colors
  )
}

# Italicize family names but leave non-taxonomic labels
# in regular font.
family_legend_labels <- function(labels) {

  as.expression(
    lapply(
      labels,
      function(label) {

        if (
          label %in% c(
            "Other",
            "Unclassified viruses"
          )
        ) {
          label
        } else {
          bquote(italic(.(label)))
        }
      }
    )
  )
}

# ============================================================
# Whole-virome composition by sample
# ============================================================

whole_sample_plot <- ggplot(
  whole_virome_plot_data,
  aes(
    x = sample_label,
    y = relative_abundance_pct,
    fill = displayed_family
  )
) +
  geom_col(
    width = 0.82
  ) +
  facet_grid(
    . ~ group_label,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = whole_virome_colors,
    breaks = family_order,
    labels = family_legend_labels
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  labs(
    x = "Sample",
    y = "Relative abundance (%)",
    fill = "Family"
  ) +
  theme_classic(base_size = 15) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 15,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 12
    ),
    axis.text.y = element_text(
      size = 12
    ),
    axis.title = element_text(
      size = 16
    ),
    legend.title = element_text(
      size = 15
    ),
    legend.text = element_text(
      size = 12
    )
  )

ggsave(
  file.path(
    output_dir,
    "whole_virome_composition_by_sample.png"
  ),
  whole_sample_plot,
  width = 14,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "whole_virome_composition_by_sample.pdf"
  ),
  whole_sample_plot,
  width = 14,
  height = 7
)

# ============================================================
# Mean composition by group
# ============================================================

group_mean_composition <- whole_virome_plot_data |>
  group_by(
    group_label,
    displayed_family
  ) |>
  summarise(
    mean_relative_abundance_pct = mean(
      relative_abundance_pct
    ),
    .groups = "drop"
  )

write_tsv(
  group_mean_composition |>
    mutate(
      group_label = as.character(group_label),
      displayed_family = as.character(
        displayed_family
      )
    ),
  file.path(
    output_dir,
    "whole_virome_group_mean_composition.tsv"
  )
)

whole_group_plot <- ggplot(
  group_mean_composition,
  aes(
    x = group_label,
    y = mean_relative_abundance_pct,
    fill = displayed_family
  )
) +
  geom_col(
    width = 0.7
  ) +
  scale_fill_manual(
    values = whole_virome_colors,
    breaks = family_order,
    labels = family_legend_labels
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  labs(
    x = "Group",
    y = "Mean relative abundance (%)",
    fill = "Family"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(
      size = 14
    ),
    axis.text.y = element_text(
      size = 12
    ),
    axis.title = element_text(
      size = 16
    ),
    legend.title = element_text(
      size = 15
    ),
    legend.text = element_text(
      size = 12
    )
  )

ggsave(
  file.path(
    output_dir,
    "whole_virome_mean_composition_by_group.png"
  ),
  whole_group_plot,
  width = 11,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "whole_virome_mean_composition_by_group.pdf"
  ),
  whole_group_plot,
  width = 11,
  height = 7
)

# ============================================================
# Descriptive bacteriophage composition
# ============================================================

categorized_counts <- read_tsv(
  file.path(
    analysis_dir,
    "viral_family_counts_with_categories.tsv"
  ),
  show_col_types = FALSE
)

phage_sample_summary <- read_tsv(
  file.path(
    analysis_dir,
    "phage_sample_summary.tsv"
  ),
  show_col_types = FALSE
)

eligible_phage_samples <- phage_sample_summary |>
  filter(
    phage_reads >= minimum_phage_reads
  ) |>
  pull(sample_id)

phage_plot_data <- categorized_counts |>
  filter(
    category == "bacteriophage",
    sample_id %in% eligible_phage_samples
  ) |>
  group_by(sample_id) |>
  mutate(
    phage_reads = sum(count),
    relative_abundance_pct = 100 * count / phage_reads
  ) |>
  ungroup() |>
  left_join(
    metadata |>
      select(
        sample_id,
        sample_label,
        group_label
      ),
    by = "sample_id"
  )

# Show the eight most abundant bacteriophage families.
top_phage_families <- phage_plot_data |>
  group_by(family) |>
  summarise(
    total_reads = sum(count),
    .groups = "drop"
  ) |>
  slice_max(
    order_by = total_reads,
    n = 8,
    with_ties = FALSE
  ) |>
  pull(family)

phage_plot_data <- phage_plot_data |>
  mutate(
    displayed_family = if_else(
      family %in% top_phage_families,
      family,
      "Other"
    )
  ) |>
  group_by(
    sample_id,
    sample_label,
    group_label,
    phage_reads,
    displayed_family
  ) |>
  summarise(
    relative_abundance_pct = sum(
      relative_abundance_pct
    ),
    .groups = "drop"
  ) |>
  mutate(
    sample_with_depth = paste0(
      sample_label,
      "\n(n=",
      comma(phage_reads),
      ")"
    )
  )

phage_family_order <- phage_plot_data |>
  group_by(displayed_family) |>
  summarise(
    total_abundance = sum(relative_abundance_pct),
    .groups = "drop"
  ) |>
  arrange(desc(total_abundance)) |>
  pull(displayed_family)

phage_family_order <- c(
  "Other",
  setdiff(
    phage_family_order,
    "Other"
  )
)

phage_plot_data <- phage_plot_data |>
  mutate(
    displayed_family = factor(
      displayed_family,
      levels = phage_family_order
    )
  )

write_tsv(
  phage_plot_data |>
    mutate(
      group_label = as.character(group_label),
      displayed_family = as.character(
        displayed_family
      )
    ),
  file.path(
    output_dir,
    "phage_composition_minimum_100_reads.tsv"
  )
)

phage_colors <- c(
  "Autotranscriptaviridae" = "#4E7189",
  "Fiersviridae" = "#4D4DFF",
  "Andersonviridae" = "#D33C32",
  "Chimalliviridae" = "#7BA05B",
  "Herelleviridae" = "#F1E487",
  "Peduoviridae" = "#C16B38",
  "Straboviridae" = "#59A9D3",
  "Casjensviridae" = "#9B59B6",
  "Other" = "#908A99"
)

missing_phage_colors <- setdiff(
  phage_family_order,
  names(phage_colors)
)

if (length(missing_phage_colors) > 0) {

  additional_phage_colors <- hue_pal()(
    length(missing_phage_colors)
  )

  names(additional_phage_colors) <-
    missing_phage_colors

  phage_colors <- c(
    phage_colors,
    additional_phage_colors
  )
}

phage_plot <- ggplot(
  phage_plot_data,
  aes(
    x = sample_with_depth,
    y = relative_abundance_pct,
    fill = displayed_family
  )
) +
  geom_col(
    width = 0.82
  ) +
  facet_grid(
    . ~ group_label,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = phage_colors,
    breaks = phage_family_order,
    labels = family_legend_labels
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  labs(
    x = "Sample and number of phage reads",
    y = "Relative abundance within bacteriophages (%)",
    fill = "Family"
  ) +
  theme_classic(base_size = 15) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 15,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 11
    ),
    axis.text.y = element_text(
      size = 12
    ),
    axis.title = element_text(
      size = 16
    ),
    legend.title = element_text(
      size = 15
    ),
    legend.text = element_text(
      size = 12
    )
  )

ggsave(
  file.path(
    output_dir,
    "phage_composition_descriptive.png"
  ),
  phage_plot,
  width = 13,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "phage_composition_descriptive.pdf"
  ),
  phage_plot,
  width = 13,
  height = 7
)

cat("\nVirome composition figures were created.\n")
cat(
  "Minimum phage depth used for the descriptive plot:",
  minimum_phage_reads,
  "\n"
)
cat(
  "Samples retained in the phage plot:",
  length(eligible_phage_samples),
  "\n"
)
cat("Results:", output_dir, "\n")
