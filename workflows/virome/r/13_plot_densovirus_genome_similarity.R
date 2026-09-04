# Plot BLASTn-based whole-genome similarity among densoviruses.

library(dplyr)
library(grid)
library(gridExtra)
library(pheatmap)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript 13_plot_densovirus_genome_similarity.R ",
    "<pair_summary.tsv> <tree_metadata.tsv> <output_directory>"
  )
}

pair_file <- args[1]
metadata_file <- args[2]
output_directory <- args[3]

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

pairs <- read_tsv(
  pair_file,
  show_col_types = FALSE
)

metadata <- read_tsv(
  metadata_file,
  show_col_types = FALSE
) %>%
  transmute(
    analysis_id = if_else(
      Name == "NC_031450.1",
      "RefSeq_NC_031450.1",
      Name
    ),
    label = case_when(
      Name == "NC_031450.1" ~ "NC_031450.1 PmDNV-JL",
      TRUE ~ sub(
        "^Densovirinae sp\\. isolate ",
        "",
        Full.Name
      )
    ),
    study = case_when(
      Dataset == "previous_study" ~ "Popov et al., 2025",
      Dataset == "current_study" &
        Group == "ABH" ~ "This study: ABH",
      Dataset == "current_study" &
        Group == "H" ~ "This study: H",
      Dataset == "current_study" &
        Group == "AAH" ~ "This study: AAH",
      Name == "NC_031450.1" ~ "Yang et al., 2016"
    )
  ) %>%
  filter(analysis_id %in% pairs$query_id)

sequence_ids <- metadata$analysis_id

# Collapse the two directional BLASTn comparisons
# into one result for each pair of genomes.
pair_summary <- pairs %>%
  filter(query_id != reference_id) %>%
  mutate(
    sequence_1 = pmin(query_id, reference_id),
    sequence_2 = pmax(query_id, reference_id)
  ) %>%
  group_by(sequence_1, sequence_2) %>%
  summarise(
    nucleotide_identity_pct =
      mean(weighted_hsp_identity_pct),
    minimum_reciprocal_coverage_pct = min(
      query_coverage_pct,
      reference_coverage_pct
    ),
    .groups = "drop"
  )

make_matrix <- function(value_column) {
  result <- matrix(
    NA_real_,
    nrow = length(sequence_ids),
    ncol = length(sequence_ids),
    dimnames = list(sequence_ids, sequence_ids)
  )

  diag(result) <- 100

  for (i in seq_len(nrow(pair_summary))) {
    id_1 <- pair_summary$sequence_1[i]
    id_2 <- pair_summary$sequence_2[i]
    value <- pair_summary[[value_column]][i]

    result[id_1, id_2] <- value
    result[id_2, id_1] <- value
  }

  result
}

identity_matrix <- make_matrix(
  "nucleotide_identity_pct"
)

coverage_matrix <- make_matrix(
  "minimum_reciprocal_coverage_pct"
)

# Save numerical results.
write_tsv(
  pair_summary,
  file.path(
    output_directory,
    "densovirus_pairwise_similarity.tsv"
  )
)

write_matrix <- function(x, filename) {
  output_table <- data.frame(
    sequence_id = rownames(x),
    x,
    check.names = FALSE
  )

  write_tsv(
    output_table,
    file.path(output_directory, filename)
  )
}

write_matrix(
  identity_matrix,
  "densovirus_whole_genome_identity_matrix.tsv"
)

write_matrix(
  coverage_matrix,
  "densovirus_minimum_reciprocal_coverage_matrix.tsv"
)

# Replace technical sequence IDs with isolate names.
display_labels <- setNames(
  metadata$label,
  metadata$analysis_id
)

rownames(identity_matrix) <-
  display_labels[rownames(identity_matrix)]

colnames(identity_matrix) <-
  display_labels[colnames(identity_matrix)]

rownames(coverage_matrix) <-
  display_labels[rownames(coverage_matrix)]

colnames(coverage_matrix) <-
  display_labels[colnames(coverage_matrix)]

annotation <- data.frame(
  Study = factor(
    metadata$study,
    levels = c(
      "Yang et al., 2016",
      "Popov et al., 2025",
      "This study: ABH",
      "This study: H",
      "This study: AAH"
    )
  ),
  row.names = metadata$label
)

annotation_colours <- list(
  Study = c(
    "Yang et al., 2016" = "#555555",
    "Popov et al., 2025" = "#E69F00",
    "This study: ABH" = "#0000FF",
    "This study: H" = "#FF0000",
    "This study: AAH" = "#BDBDBD"
  )
)

# Use nucleotide identity for a common row and column order.
clustering_tree <- hclust(
  as.dist(100 - identity_matrix),
  method = "average"
)

plot_heatmap <- function(
  data_matrix,
  legend_title,
  colours,
  breaks,
  legend_breaks,
  output_stem
) {
  numbers <- matrix(
    sprintf("%.1f", data_matrix),
    nrow = nrow(data_matrix),
    dimnames = dimnames(data_matrix)
  )

  common_arguments <- list(
      mat = data_matrix,
      color = colours,
      breaks = breaks,
      cluster_rows = clustering_tree,
      cluster_cols = clustering_tree,
      annotation_row = annotation,
      annotation_col = annotation,
      annotation_colors = annotation_colours,
      annotation_names_row = FALSE,
      annotation_names_col = FALSE,
      display_numbers = numbers,
      fontsize = 14,
      fontsize_row = 11,
      fontsize_col = 11,
      fontsize_number = 7.5,
      angle_col = 45,
      border_color = "white",
      legend_breaks = legend_breaks,
      legend_labels = sprintf("%.1f", legend_breaks),
      treeheight_row = 65,
      treeheight_col = 65,
      silent = TRUE
  )

  legend_arguments <- common_arguments
  legend_arguments$fontsize <- 16

  heatmap_with_legends <- do.call(
    pheatmap,
    c(
      legend_arguments,
      list(
        legend = TRUE,
        annotation_legend = TRUE
      )
    )
  )

  heatmap_without_legends <- do.call(
    pheatmap,
    c(
      common_arguments,
      list(
        legend = FALSE,
        annotation_legend = FALSE
      )
    )
  )

  extract_grob <- function(gtable_object, grob_name) {
    position <- which(
      gtable_object$layout$name == grob_name
    )

    if (length(position) == 0) {
      return(nullGrob())
    }

    gtable_object$grobs[[position[1]]]
  }

  colour_legend <- extract_grob(
    heatmap_with_legends$gtable,
    "legend"
  )

  study_legend <- extract_grob(
    heatmap_with_legends$gtable,
    "annotation_legend"
  )

  title_height_mm <- if (
    grepl("\n", legend_title, fixed = TRUE)
  ) {
    12
  } else {
    7
  }

  colour_legend_height_mm <-
    title_height_mm + 4 + 27

  legend_block_height_mm <-
    colour_legend_height_mm + 42

  titled_colour_legend <- arrangeGrob(
    textGrob(
      legend_title,
      x = unit(0, "npc"),
      just = "left",
      gp = gpar(
        fontsize = 16,
        fontface = "bold",
        lineheight = 0.9
      )
    ),
    nullGrob(),
    colour_legend,
    ncol = 1,
    heights = unit(
      c(title_height_mm, 4, 27),
      "mm"
    )
  )

  legend_block <- arrangeGrob(
    titled_colour_legend,
    study_legend,
    ncol = 1,
    heights = unit(
      c(colour_legend_height_mm, 42),
      "mm"
    )
  )

  centred_legend_panel <- arrangeGrob(
    nullGrob(),
    legend_block,
    nullGrob(),
    ncol = 1,
    heights = unit.c(
      unit(1, "null"),
      unit(legend_block_height_mm + 4, "mm"),
      unit(1, "null")
    )
  )

  complete_plot <- arrangeGrob(
    heatmap_without_legends$gtable,
    centred_legend_panel,
    ncol = 2,
    widths = unit.c(
      unit(1, "null"),
      unit(62, "mm")
    )
  )

  for (extension in c("png", "pdf")) {
    output_file <- file.path(
      output_directory,
      paste0(output_stem, ".", extension)
    )

    if (extension == "png") {
      png(
        output_file,
        width = 15,
        height = 12.5,
        units = "in",
        res = 300
      )
    } else {
      pdf(
        output_file,
        width = 15,
        height = 12.5,
        useDingbats = FALSE
      )
    }

    grid.newpage()
    grid.draw(complete_plot)
    dev.off()
  }
}

identity_minimum <-
  floor(min(identity_matrix) * 2) / 2

plot_heatmap(
  identity_matrix,
  "Identity (%)",
  colorRampPalette(
    c("#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8")
  )(100),
  seq(
    identity_minimum,
    100,
    length.out = 101
  ),
  seq(
    ceiling(identity_minimum),
    100,
    by = 1
  ),
  "densovirus_whole_genome_identity_heatmap"
)

coverage_minimum <- floor(min(coverage_matrix))

plot_heatmap(
  coverage_matrix,
  "Reciprocal\ncoverage (%)",
  colorRampPalette(
    c("#F7FBFF", "#9ECAE1", "#08519C")
  )(100),
  seq(
    coverage_minimum,
    100,
    length.out = 101
  ),
  unique(
    c(
      coverage_minimum,
      90,
      95,
      100
    )
  ),
  "densovirus_minimum_reciprocal_coverage_heatmap"
)

message("Densovirus similarity heatmaps completed.")
message("Sequences: ", length(sequence_ids))
message(
  "Unique non-self pairs: ",
  nrow(pair_summary)
)
message(
  "Minimum nucleotide identity: ",
  round(min(identity_matrix), 3),
  "%"
)
message(
  "Minimum reciprocal coverage: ",
  round(min(coverage_matrix), 3),
  "%"
)
message("Results: ", output_directory)
