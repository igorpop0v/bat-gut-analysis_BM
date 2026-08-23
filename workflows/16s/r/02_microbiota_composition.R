# Summarize and visualize microbiota composition at five taxonomic ranks.

library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(patchwork)

# Local input and output paths.
source("workflows/16s/r/local_paths.R")

# Input files.
counts_file <- file.path(
  input_dir,
  "asv_counts_unrarefied.tsv"
)

metadata_file <- file.path(
  input_dir,
  "sample_metadata.tsv"
)

taxonomy_file <- file.path(
  results_dir,
  "taxonomy",
  "asv_taxonomy_parsed_to_genus.tsv"
)

# Output directory.
output_dir <- file.path(
  results_dir,
  "microbiota_composition"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Order and labels of the experimental groups.
group_levels <- c(
  "active",
  "hibernation",
  "post_hibernation"
)

# ABH: active before hibernation
# H: hibernation
# AAH: active after hibernation
group_labels <- c(
  active = "ABH",
  hibernation = "H",
  post_hibernation = "AAH"
)

# Taxa below these thresholds in every group are combined as "Other".
# Values are proportions: 0.01 = 1%, 0.05 = 5%.
rank_thresholds <- c(
  phylum = 0.01,
  class = 0.01,
  order = 0.01,
  family = 0.05,
  genus = 0.01
)

# Read the unrarefied ASV count table.
counts <- read_tsv(
  counts_file,
  skip = 1,
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  rename(feature_id = 1)

# Read sample metadata.
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

# Read parsed GTDB taxonomy.
taxonomy <- read_tsv(
  taxonomy_file,
  show_col_types = FALSE
) |>
  mutate(
    # Add the family name to GTDB placeholder genera.
    genus = case_when(
      grepl("^(JALENY|GCF-)", genus) ~ paste(family, genus),
      TRUE ~ genus
    )
  )

# Convert the ASV table to long format and calculate
# relative abundance separately within each sample.
counts_long <- counts |>
  pivot_longer(
    -feature_id,
    names_to = "sample_id",
    values_to = "count"
  ) |>
  left_join(
    metadata,
    by = "sample_id"
  ) |>
  left_join(
    taxonomy,
    by = "feature_id"
  ) |>
  group_by(sample_id) |>
  mutate(
    relative_abundance = count / sum(count)
  ) |>
  ungroup()

# Summarize one taxonomic rank.
summarize_rank <- function(rank_name, threshold) {

  rank_summary <- counts_long |>
    transmute(
      sample_id,
      group,
      taxon = .data[[rank_name]],
      relative_abundance
    ) |>
    group_by(
      sample_id,
      group,
      taxon
    ) |>
    summarise(
      relative_abundance = sum(relative_abundance),
      .groups = "drop"
    ) |>
    group_by(
      group,
      taxon
    ) |>
    summarise(
      mean_relative_abundance = mean(relative_abundance),
      .groups = "drop"
    )

  # Identify taxa that do not reach the threshold in any group.
  pooled_taxa <- rank_summary |>
    group_by(taxon) |>
    summarise(
      maximum_group_mean = max(mean_relative_abundance),
      .groups = "drop"
    ) |>
    filter(maximum_group_mean < threshold) |>
    pull(taxon)

  # Combine low-abundance taxa as "Other".
  rank_summary |>
    mutate(
      taxon = if_else(
        taxon %in% pooled_taxa,
        "Other",
        taxon
      )
    ) |>
    group_by(
      group,
      taxon
    ) |>
    summarise(
      mean_relative_abundance_percent =
        100 * sum(mean_relative_abundance),
      .groups = "drop"
    ) |>
    mutate(rank = rank_name)
}

# Process all five ranks.
rank_summaries <- imap(
  rank_thresholds,
  ~ summarize_rank(.y, .x)
)

# Save the underlying numerical data.
composition_table <- bind_rows(rank_summaries) |>
  select(
    rank,
    group,
    taxon,
    mean_relative_abundance_percent
  )

write_tsv(
  composition_table,
  file.path(
    output_dir,
    "mean_relative_abundance_by_group.tsv"
  )
)

# Prepare taxon labels for the figure.
taxon_labels <- function(labels, rank_name) {

  lapply(labels, function(label) {

    if (label == "Other") {
      return(label)
    }

    # Replace underscores only in displayed labels.
    display_label <- gsub(
      "_",
      " ",
      label,
      fixed = TRUE
    )

    # Italicize the family but not the GTDB placeholder code.
    if (
      rank_name == "genus" &&
      grepl(
        " (JALENY|GCF-)[^ ]+$",
        display_label
      )
    ) {

      family_name <- sub(
        " ([^ ]+)$",
        "",
        display_label
      )

      genus_code <- sub(
        "^.* ",
        "",
        display_label
      )

      return(
        bquote(
          italic(.(family_name)) ~ .(genus_code)
        )
      )
    }

    # Keep "Unclassified" in regular font and italicize
    # only the taxonomic name that follows it.
    if (
      rank_name %in% c("family", "genus") &&
      startsWith(display_label, "Unclassified ")
    ) {
      
      parent_taxon <- sub(
        "^Unclassified ",
        "",
        display_label
      )
      
      return(
        bquote(
          "Unclassified " * italic(.(parent_taxon))
        )
      )
    }
    
    # Italicize classified family- and genus-level labels.
    if (rank_name %in% c("family", "genus")) {
      return(
        bquote(
          italic(.(display_label))
        )
      )
    }

    display_label
  })
}

# Create a heatmap for one taxonomic rank.
plot_rank_heatmap <- function(
    rank_data,
    rank_name,
    fill_max,
    taxa_on_x = FALSE) {

  # Order taxa by their mean abundance across groups.
  taxon_order <- rank_data |>
    group_by(taxon) |>
    summarise(
      overall_mean = mean(
        mean_relative_abundance_percent
      ),
      .groups = "drop"
    ) |>
    arrange(desc(overall_mean)) |>
    pull(taxon)

  # Always place "Other" last.
  taxon_order <- c(
    setdiff(taxon_order, "Other"),
    intersect("Other", taxon_order)
  )

  # Vertical heatmaps for phylum, class, order and family.
  if (!taxa_on_x) {

    plot_data <- rank_data |>
      mutate(
        taxon = factor(
          taxon,
          levels = rev(taxon_order)
        )
      )

    return(
      ggplot(
        plot_data,
        aes(
          x = group,
          y = taxon,
          fill = mean_relative_abundance_percent
        )
      ) +
        geom_tile(
          color = "white",
          linewidth = 0.4
        ) +
        geom_text(
          aes(
            label = sprintf(
              "%.1f",
              mean_relative_abundance_percent
            )
          ),
          size = 3.2
        ) +
        scale_x_discrete(
          labels = group_labels,
          drop = FALSE
        ) +
        scale_y_discrete(
          labels = function(x) {
            taxon_labels(x, rank_name)
          }
        ) +
        scale_fill_gradient(
          name = "Mean RA (%)",
          low = "#FFFFFF",
          high = "#D9A45B",
          limits = c(0, fill_max),
          oob = scales::squish,
          guide = guide_colorbar(
            barheight = grid::unit(60, "mm"),
            barwidth = grid::unit(7, "mm")
          )
        ) +
        labs(
          x = NULL,
          y = NULL
        ) +
        theme_minimal(base_size = 10) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(
            color = "black",
            size = 12,
            lineheight = 0.9
          ),
          axis.text.y = element_text(
            color = "black"
          ),
          legend.title = element_text(
            hjust = 0.5
          ),
          plot.margin = margin(
            5.5,
            8,
            5.5,
            5.5
          )
        )
    )
  }

  # Wide horizontal heatmap for genus.
  plot_data <- rank_data |>
    mutate(
      taxon = factor(
        taxon,
        levels = taxon_order
      )
    )

  ggplot(
    plot_data,
    aes(
      x = taxon,
      y = group,
      fill = mean_relative_abundance_percent
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.4
    ) +
    geom_text(
      aes(
        label = sprintf(
          "%.1f",
          mean_relative_abundance_percent
        )
      ),
      size = 3.2
    ) +
    scale_x_discrete(
      labels = function(x) {
        taxon_labels(x, rank_name)
      }
    ) +
    scale_y_discrete(
      labels = group_labels,
      drop = FALSE
    ) +
    scale_fill_gradient(
      name = "Mean RA (%)",
      low = "#FFFFFF",
      high = "#D9A45B",
      limits = c(0, fill_max),
      oob = scales::squish,
      guide = guide_colorbar(
        barheight = grid::unit(60, "mm"),
        barwidth = grid::unit(7, "mm")
      )
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 35,
        hjust = 1,
        vjust = 1,
        color = "black"
      ),
      axis.text.y = element_text(
        color = "black",
        size = 12,
        lineheight = 0.9
      ),
      legend.title = element_text(
        hjust = 0.5
      ),
      plot.margin = margin(
        5.5,
        8,
        5.5,
        5.5
      )
    )
}

# Use the same colour scale in all panels.
fill_max <- composition_table |>
  summarise(
    value = 10 * ceiling(
      max(mean_relative_abundance_percent) / 10
    )
  ) |>
  pull(value)

# Create individual panels.
phylum_plot <- plot_rank_heatmap(
  rank_summaries$phylum,
  "phylum",
  fill_max
)

class_plot <- plot_rank_heatmap(
  rank_summaries$class,
  "class",
  fill_max
)

order_plot <- plot_rank_heatmap(
  rank_summaries$order,
  "order",
  fill_max
)

family_plot <- plot_rank_heatmap(
  rank_summaries$family,
  "family",
  fill_max
)

genus_plot <- plot_rank_heatmap(
  rank_summaries$genus,
  "genus",
  fill_max,
  taxa_on_x = TRUE
)

# Combine the five panels.
combined_plot <- (
  (phylum_plot | class_plot) /
    (order_plot | family_plot) /
    genus_plot
) +
  plot_layout(
    heights = c(1, 1.2, 1.25),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "right",
    legend.title = element_text(
      size = 12,
      hjust = 0.5
    ),
    legend.text = element_text(
      size = 11
    )
  )

# Save raster and vector versions.
ggsave(
  file.path(
    output_dir,
    "microbiota_composition_heatmap.png"
  ),
  combined_plot,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "microbiota_composition_heatmap.pdf"
  ),
  combined_plot,
  width = 14,
  height = 10,
  units = "in",
  bg = "white"
)

message("Microbiota-composition analysis finished.")
message("Output directory: ", output_dir)
