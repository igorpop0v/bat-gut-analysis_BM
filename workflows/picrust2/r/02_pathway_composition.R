# Composition of predicted MetaCyc pathways in the three groups.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

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
  "pathway_composition"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Number of pathways displayed in the heatmap.
top_n <- 25

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

# Read MetaCyc pathway abundances.
pathway_table <- read_tsv(
  pathway_file,
  show_col_types = FALSE,
  name_repair = "minimal"
)

# Read biological sample metadata.
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

sample_columns <- intersect(
  names(pathway_table),
  metadata$sample_id
)

if (length(sample_columns) != nrow(metadata)) {
  stop(
    "The pathway table and metadata do not contain ",
    "the same biological samples."
  )
}

# Convert the pathway table into long format.
pathway_long <- pathway_table |>
  mutate(
    description = if_else(
      is.na(description) |
        description == "not_found",
      pathway,
      description
    )
  ) |>
  pivot_longer(
    cols = all_of(sample_columns),
    names_to = "sample_id",
    values_to = "predicted_abundance"
  ) |>
  left_join(
    metadata,
    by = "sample_id"
  )

# Convert predicted pathway abundances to relative abundance
# separately within every sample.
pathway_relative <- pathway_long |>
  group_by(sample_id) |>
  mutate(
    total_predicted_abundance =
      sum(predicted_abundance),
    relative_abundance =
      predicted_abundance /
      total_predicted_abundance
  ) |>
  ungroup()

if (any(
  !is.finite(pathway_relative$relative_abundance)
)) {
  stop(
    "Non-finite relative pathway abundances were produced."
  )
}

# Save the complete sample-level relative abundance table.
write_tsv(
  pathway_relative |>
    select(
      sample_id,
      group,
      pathway,
      description,
      predicted_abundance,
      relative_abundance
    ) |>
    mutate(
      group = as.character(group)
    ),
  file.path(
    output_dir,
    "pathway_relative_abundance_by_sample.tsv"
  )
)

# Calculate mean relative abundance within each group.
group_means <- pathway_relative |>
  group_by(
    pathway,
    description,
    group
  ) |>
  summarise(
    mean_relative_abundance =
      mean(relative_abundance),
    .groups = "drop"
  ) |>
  mutate(
    mean_relative_abundance_percent =
      100 * mean_relative_abundance
  )

write_tsv(
  group_means |>
    mutate(
      group = as.character(group)
    ),
  file.path(
    output_dir,
    "pathway_mean_relative_abundance_by_group.tsv"
  )
)

# Calculate overall means across all samples.
overall_means <- pathway_relative |>
  group_by(
    pathway,
    description
  ) |>
  summarise(
    overall_mean_relative_abundance =
      mean(relative_abundance),
    .groups = "drop"
  )

# Rank pathways according to their maximum group mean.
pathway_ranking <- group_means |>
  group_by(
    pathway,
    description
  ) |>
  summarise(
    maximum_group_mean =
      max(mean_relative_abundance),
    .groups = "drop"
  ) |>
  left_join(
    overall_means,
    by = c(
      "pathway",
      "description"
    )
  ) |>
  arrange(
    desc(maximum_group_mean),
    desc(overall_mean_relative_abundance)
  )

selected_pathways <- pathway_ranking |>
  slice_head(n = top_n)

write_tsv(
  selected_pathways,
  file.path(
    output_dir,
    "pathways_selected_for_heatmap.tsv"
  )
)

# Prepare the selected pathways for clustering.
heatmap_data <- group_means |>
  semi_join(
    selected_pathways,
    by = c(
      "pathway",
      "description"
    )
  )

cluster_table <- heatmap_data |>
  select(
    pathway,
    group,
    mean_relative_abundance_percent
  ) |>
  pivot_wider(
    names_from = group,
    values_from = mean_relative_abundance_percent
  )

cluster_matrix <- cluster_table |>
  select(
    all_of(group_levels)
  ) |>
  as.matrix()

rownames(cluster_matrix) <- cluster_table$pathway

pathway_order <- rownames(
  cluster_matrix
)[
  hclust(
    dist(cluster_matrix),
    method = "complete"
  )$order
]

# Combine the readable description with the stable MetaCyc ID.
pathway_labels <- selected_pathways |>
  mutate(
    display_label = paste0(
      description,
      " (",
      pathway,
      ")"
    )
  )

label_map <- setNames(
  pathway_labels$display_label,
  pathway_labels$pathway
)

wrap_pathway_labels <- function(pathway_ids) {
  vapply(
    pathway_ids,
    function(pathway_id) {
      paste(
        strwrap(
          label_map[[pathway_id]],
          width = 52
        ),
        collapse = "\n"
      )
    },
    character(1)
  )
}

heatmap_data <- heatmap_data |>
  mutate(
    pathway = factor(
      pathway,
      levels = rev(pathway_order)
    ),
    group = factor(
      group,
      levels = group_levels
    )
  )

fill_max <- max(
  heatmap_data$mean_relative_abundance_percent
)

pathway_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = group,
    y = pathway,
    fill = mean_relative_abundance_percent
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        mean_relative_abundance_percent
      )
    ),
    size = 4.2
  ) +
  scale_x_discrete(
    labels = group_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    labels = wrap_pathway_labels
  ) +
  scale_fill_gradient(
    name = "Mean RA (%)",
    low = "#FFFFFF",
    high = "#D9A45B",
    limits = c(
      0,
      fill_max
    ),
    oob = scales::squish,
    guide = guide_colorbar(
      barheight = grid::unit(
        65,
        "mm"
      ),
      barwidth = grid::unit(
        8,
        "mm"
      )
    )
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_size = 15
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 16,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11.5,
      color = "black"
    ),
    legend.title = element_text(
      size = 15
    ),
    legend.text = element_text(
      size = 13
    ),
    plot.margin = margin(
      10,
      20,
      10,
      10
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    paste0(
      "predicted_pathway_composition_top",
      top_n,
      ".png"
    )
  ),
  plot = pathway_heatmap,
  width = 14,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("")
message("Predicted pathway composition analysis completed.")
message("Results: ", output_dir)
message("")
print(selected_pathways)
