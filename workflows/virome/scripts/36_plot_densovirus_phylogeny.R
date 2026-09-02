#!/usr/bin/env Rscript

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

tree_file <- args[1]
manifest_file <- args[2]
custom_labels_file <- args[3]
output_dir <- args[4]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

minimum_alrt <- 80
minimum_ufboot <- 95

outgroup_ids <- c(
  "NC_011545.2",
  "NC_077018.1",
  "NC_015115.1",
  "NC_011317.1",
  "NC_012636.1"
)

tree <- read.tree(tree_file)

if (!all(outgroup_ids %in% tree$tip.label)) {
  stop("One or more outgroup sequences are absent.")
}

tree <- root(
  tree,
  outgroup = outgroup_ids,
  resolve.root = TRUE
)

tree <- ladderize(
  tree,
  right = FALSE
)

manifest <- read.delim(
  manifest_file,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

custom_labels <- read.delim(
  custom_labels_file,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

tip_metadata <- manifest %>%
  filter(analysis_id %in% tree$tip.label) %>%
  mutate(
    display_label = sub(
      ", complete.*$",
      "",
      description
    ),
    tip_class = case_when(
      dataset == "current_study" &
        group == "ABH" ~ "Current study: ABH",

      dataset == "current_study" &
        group == "H" ~ "Current study: H",

      dataset == "current_study" &
        group == "AAH" ~ "Current study: AAH",

      dataset == "previous_study" ~
        "Previous study",

      TRUE ~ "External reference"
    )
  ) %>%
  select(
    analysis_id,
    display_label,
    tip_class,
    dataset,
    group
  ) %>%
  left_join(
    custom_labels %>%
      rename(custom_display_label = display_label),
    by = "analysis_id"
  ) %>%
  mutate(
    display_label = if_else(
      !is.na(custom_display_label),
      custom_display_label,
      display_label
    )
  ) %>%
  transmute(
    label = analysis_id,
    analysis_id,
    display_label,
    tip_class,
    dataset,
    group
  )

write.table(
  tip_metadata,
  file.path(
    output_dir,
    "densovirus_tree_labels_resolved.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

parse_support <- function(values, position) {

  vapply(
    values,
    function(value) {

      if (is.na(value) || value == "") {
        return(NA_real_)
      }

      parts <- strsplit(
        value,
        "/",
        fixed = TRUE
      )[[1]]

      if (length(parts) < position) {
        return(NA_real_)
      }

      suppressWarnings(
        as.numeric(parts[position])
      )
    },
    numeric(1)
  )
}

add_support_information <- function(plot_object) {

  plot_object$data$sh_alrt <- parse_support(
    plot_object$data$label,
    1
  )

  plot_object$data$ufboot <- parse_support(
    plot_object$data$label,
    2
  )

  plot_object$data$strong_support <- (
    !plot_object$data$isTip &
    !is.na(plot_object$data$sh_alrt) &
    !is.na(plot_object$data$ufboot) &
    plot_object$data$sh_alrt >= minimum_alrt &
    plot_object$data$ufboot >= minimum_ufboot
  )

  plot_object$data$support_label <- ifelse(
    plot_object$data$strong_support,
    paste0(
      format(
        plot_object$data$sh_alrt,
        trim = TRUE
      ),
      "/",
      format(
        plot_object$data$ufboot,
        trim = TRUE
      )
    ),
    NA_character_
  )

  plot_object
}

class_order <- c(
  "External reference",
  "Previous study",
  "Current study: ABH",
  "Current study: H",
  "Current study: AAH"
)

class_colors <- c(
  "External reference" = "#4D4D4D",
  "Previous study" = "#E69F00",
  "Current study: ABH" = "#0000FF",
  "Current study: H" = "#FF0000",
  "Current study: AAH" = "#BDBDBD"
)

class_shapes <- c(
  "External reference" = 16,
  "Previous study" = 15,
  "Current study: ABH" = 16,
  "Current study: H" = 16,
  "Current study: AAH" = 16
)

study_ids <- grep(
  "^(Current|Previous)_",
  tree$tip.label,
  value = TRUE
)

focal_node <- getMRCA(
  tree,
  study_ids
)

outgroup_node <- getMRCA(
  tree,
  outgroup_ids
)

# Complete reference tree.
full_tree_plot <- ggtree(
  tree,
  linewidth = 0.35
)

full_tree_plot <- add_support_information(
  full_tree_plot
)

full_max_x <- max(
  full_tree_plot$data$x,
  na.rm = TRUE
)

full_tree_plot <- full_tree_plot %<+%
  tip_metadata

full_tree_plot <- full_tree_plot +
  geom_hilight(
    node = outgroup_node,
    fill = "#D9D9D9",
    alpha = 0.25
  ) +
  geom_hilight(
    node = focal_node,
    fill = "#FFF2CC",
    alpha = 0.35
  ) +
  geom_point2(
    aes(
      subset = !isTip & strong_support
    ),
    shape = 16,
    size = 0.9,
    colour = "black"
  ) +
  geom_tippoint(
    aes(
      colour = tip_class,
      shape = tip_class
    ),
    size = 1.8
  ) +
  geom_tiplab(
    aes(label = display_label),
    align = TRUE,
    linetype = "dotted",
    linesize = 0.2,
    size = 2.4,
    offset = 0.01
  ) +
  scale_colour_manual(
    values = class_colors,
    breaks = class_order,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = class_shapes,
    breaks = class_order,
    drop = FALSE
  ) +
  xlim(
    0,
    full_max_x * 1.85
  ) +
  labs(
    x = "Amino-acid substitutions per site",
    colour = NULL,
    shape = NULL
  ) +
  theme_tree2() +
  theme(
    text = element_text(size = 11),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(
      8,
      20,
      8,
      8
    )
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE
    ),
    shape = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )

# Focal clade containing the study sequences and closest references.
focal_tip_ids <- c(
  study_ids,
  "NC_031450.1",
  "MW628494.1"
)

focal_tip_ids <- intersect(
  focal_tip_ids,
  tree$tip.label
)

focal_tree <- keep.tip(
  tree,
  focal_tip_ids
)

focal_tree <- ladderize(
  focal_tree,
  right = FALSE
)

focal_tree_plot <- ggtree(
  focal_tree,
  linewidth = 0.5
)

focal_tree_plot <- add_support_information(
  focal_tree_plot
)

focal_max_x <- max(
  focal_tree_plot$data$x,
  na.rm = TRUE
)

focal_tree_plot <- focal_tree_plot %<+%
  tip_metadata

focal_tree_plot <- focal_tree_plot +
  geom_point2(
    aes(
      subset = !isTip & strong_support
    ),
    shape = 16,
    size = 1.4,
    colour = "black"
  ) +
  geom_text2(
    aes(
      subset = !isTip & strong_support,
      label = support_label
    ),
    size = 2.8,
    hjust = -0.15,
    vjust = -0.3
  ) +
  geom_tippoint(
    aes(
      colour = tip_class,
      shape = tip_class
    ),
    size = 2.4
  ) +
  geom_tiplab(
    aes(label = display_label),
    align = TRUE,
    linetype = "dotted",
    linesize = 0.25,
    size = 3.1,
    offset = 0.0002
  ) +
  scale_colour_manual(
    values = class_colors,
    breaks = class_order,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = class_shapes,
    breaks = class_order,
    drop = FALSE
  ) +
  xlim(
    0,
    focal_max_x * 2.25
  ) +
  labs(
    x = "Amino-acid substitutions per site",
    colour = NULL,
    shape = NULL
  ) +
  theme_tree2() +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom",
    plot.margin = margin(
      8,
      20,
      8,
      8
    )
  )

ggsave(
  file.path(
    output_dir,
    "densovirus_partitioned_full_tree.png"
  ),
  full_tree_plot,
  width = 15,
  height = 18,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "densovirus_partitioned_full_tree.pdf"
  ),
  full_tree_plot,
  width = 15,
  height = 18,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(
    output_dir,
    "densovirus_partitioned_focal_clade.png"
  ),
  focal_tree_plot,
  width = 14,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "densovirus_partitioned_focal_clade.pdf"
  ),
  focal_tree_plot,
  width = 14,
  height = 9,
  units = "in",
  device = cairo_pdf
)

cat("Densovirus phylogeny figures created.\n")
cat("Results:", output_dir, "\n")
