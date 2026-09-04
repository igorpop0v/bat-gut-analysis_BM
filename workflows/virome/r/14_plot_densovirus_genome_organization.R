#!/usr/bin/env Rscript

# Plot the organization of validated densovirus CDS.
# ORFs without reliable coordinates are reported as unresolved and are not drawn.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript 14_plot_densovirus_genome_organization.R ",
    "<densovirus_annotation_dir> <densovirus_comparative_dir> <output_dir>"
  )
}

annotation_dir <- normalizePath(args[1], mustWork = TRUE)
comparative_dir <- normalizePath(args[2], mustWork = TRUE)
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

current_gff <- file.path(
  annotation_dir,
  "provisional_annotations",
  "current_study_densovirus_provisional_annotations.gff3"
)

previous_gff <- file.path(
  annotation_dir,
  "previous_study_provisional_annotations",
  "previous_study_densovirus_provisional_annotations.gff3"
)

current_manifest_file <- file.path(
  annotation_dir,
  "curated_annotation_set",
  "current_study_annotation_manifest.tsv"
)

previous_manifest_file <- file.path(
  annotation_dir,
  "previous_study_annotation_set",
  "previous_study_annotation_manifest.tsv"
)

current_unresolved_file <- file.path(
  annotation_dir,
  "curated_annotation_set",
  "current_study_unresolved_orfs.tsv"
)

previous_unresolved_file <- file.path(
  annotation_dir,
  "previous_study_annotation_set",
  "previous_study_unresolved_orfs.tsv"
)

refseq_manifest_file <- file.path(
  comparative_dir,
  "refseq_annotation",
  "NC_031450.1_protein_manifest.tsv"
)

read_gff_cds <- function(path) {
  read_tsv(
    path,
    comment = "#",
    col_names = c(
      "analysis_id", "source", "type", "start", "end",
      "score", "strand", "phase", "attributes"
    ),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) |>
    filter(type == "CDS") |>
    transmute(
      analysis_id,
      start = as.numeric(start),
      end = as.numeric(end),
      strand,
      orf = str_match(attributes, "(?:^|;)Name=([^;]+)")[, 2]
    )
}

parse_refseq_location <- function(location) {
  coordinates <- str_match(location, "([0-9]+)\\.\\.([0-9]+)")

  tibble(
    start = as.numeric(coordinates[, 2]),
    end = as.numeric(coordinates[, 3]),
    strand = if_else(str_detect(location, "^complement\\("), "-", "+")
  )
}

current_manifest <- read_tsv(
  current_manifest_file,
  show_col_types = FALSE
) |>
  transmute(
    analysis_id,
    isolate,
    group,
    length_nt,
    study = paste0("This study: ", group)
  )

previous_manifest <- read_tsv(
  previous_manifest_file,
  show_col_types = FALSE
) |>
  transmute(
    analysis_id,
    isolate,
    group = NA_character_,
    length_nt,
    study = "Popov et al., 2025"
  )

refseq_manifest <- read_tsv(
  refseq_manifest_file,
  show_col_types = FALSE
)

refseq_features <- bind_cols(
  refseq_manifest |>
    transmute(
      analysis_id = "RefSeq_NC_031450.1",
      orf = product
    ),
  parse_refseq_location(refseq_manifest$location)
) |>
  select(analysis_id, start, end, strand, orf)

refseq_sequence <- tibble(
  analysis_id = "RefSeq_NC_031450.1",
  isolate = "NC_031450.1 PmDNV-JL",
  group = NA_character_,
  length_nt = 5166,
  study = "Yang et al., 2016"
)

features <- bind_rows(
  refseq_features,
  read_gff_cds(previous_gff),
  read_gff_cds(current_gff)
)

sequences <- bind_rows(
  refseq_sequence,
  previous_manifest,
  current_manifest
)

sequence_order <- c(
  "RefSeq_NC_031450.1",
  paste0("Previous_D", 1:5),
  paste0("Previous_P", 1:5),
  "Current_Nr11",
  "Current_Nr21", "Current_Nr22",
  "Current_Nr6", "Current_Nr7", "Current_Nr8", "Current_Nr9"
)

sequences <- sequences |>
  mutate(
    analysis_id = factor(analysis_id, levels = sequence_order),
    y = rev(seq_len(n()))
  ) |>
  arrange(analysis_id)

unresolved <- bind_rows(
  read_tsv(current_unresolved_file, show_col_types = FALSE),
  read_tsv(previous_unresolved_file, show_col_types = FALSE)
) |>
  group_by(analysis_id) |>
  summarise(
    unresolved_orfs = paste(sort(unique(orf)), collapse = ", "),
    .groups = "drop"
  )

sequences <- sequences |>
  left_join(unresolved, by = "analysis_id") |>
  mutate(
    unresolved_label = if_else(
      is.na(unresolved_orfs),
      "",
      paste0("Unresolved: ", unresolved_orfs)
    )
  )

lane_offset <- c(
  ORF5 = 0.18,
  ORF3 = 0.18,
  ORF4 = 0.38,
  ORF1 = -0.18,
  ORF2 = -0.38
)

features <- features |>
  left_join(
    sequences |> select(analysis_id, isolate, study, y),
    by = "analysis_id"
  ) |>
  mutate(
    feature_id = paste(analysis_id, orf, sep = "__"),
    centre_y = y + unname(lane_offset[orf])
  )

make_arrow <- function(feature_row, head_length = 115) {
  start <- feature_row$start
  end <- feature_row$end
  centre_y <- feature_row$centre_y
  strand <- feature_row$strand
  body_half <- 0.075
  head_half <- 0.155

  if (strand == "+") {
    head_base <- max(start, end - head_length)
    x <- c(start, head_base, head_base, end, head_base, head_base, start)
  } else {
    head_base <- min(end, start + head_length)
    x <- c(end, head_base, head_base, start, head_base, head_base, end)
  }

  y <- c(
    centre_y - body_half,
    centre_y - body_half,
    centre_y - head_half,
    centre_y,
    centre_y + head_half,
    centre_y + body_half,
    centre_y + body_half
  )

  tibble(
    feature_id = feature_row$feature_id,
    analysis_id = feature_row$analysis_id,
    orf = feature_row$orf,
    vertex = seq_along(x),
    x = x,
    y = y
  )
}

arrow_polygons <- bind_rows(
  lapply(seq_len(nrow(features)), function(index) {
    make_arrow(features[index, ])
  })
)

orf_colours <- c(
  ORF5 = "#E69F00",
  ORF3 = "#0072B2",
  ORF4 = "#56B4E9",
  ORF1 = "#009E73",
  ORF2 = "#CC79A7"
)

study_colours <- c(
  "Yang et al., 2016" = "#595959",
  "Popov et al., 2025" = "#E69F00",
  "This study: ABH" = "#0000FF",
  "This study: H" = "#FF0000",
  "This study: AAH" = "#BDBDBD"
)

plot <- ggplot() +
  geom_hline(
    yintercept = c(17.5, 7.5, 6.5, 4.5),
    colour = "#D9D9D9",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_segment(
    data = sequences,
    aes(x = 0, xend = length_nt, y = y, yend = y),
    colour = "#404040",
    linewidth = 0.6
  ) +
  geom_polygon(
    data = arrow_polygons,
    aes(x = x, y = y, group = feature_id, fill = orf),
    colour = "black",
    linewidth = 0.25
  ) +
  geom_point(
    data = sequences,
    aes(x = -210, y = y, colour = study),
    shape = 15,
    size = 4.6
  ) +
  geom_text(
    data = sequences |> filter(unresolved_label != ""),
    aes(x = 5850, y = y, label = unresolved_label),
    hjust = 0,
    colour = "#595959",
    size = 3.8
  ) +
  scale_fill_manual(
    name = "Validated CDS",
    values = orf_colours,
    breaks = names(orf_colours)
  ) +
  scale_colour_manual(
    name = "Study",
    values = study_colours,
    breaks = names(study_colours)
  ) +
  scale_x_continuous(
    name = "Genome position (kb)",
    breaks = seq(0, 6000, by = 1000),
    labels = seq(0, 6, by = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = sequences$y,
    labels = sequences$isolate,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  coord_cartesian(xlim = c(-320, 6750), clip = "off") +
  guides(
    fill = guide_legend(order = 1, override.aes = list(colour = "black")),
    colour = guide_legend(order = 2, override.aes = list(shape = 15, size = 5))
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13),
    legend.key.height = unit(0.65, "cm"),
    plot.margin = margin(12, 20, 12, 12)
  )

ggsave(
  file.path(output_dir, "densovirus_genome_organization.pdf"),
  plot,
  width = 15,
  height = 12.5,
  units = "in"
)

ggsave(
  file.path(output_dir, "densovirus_genome_organization.png"),
  plot,
  width = 15,
  height = 12.5,
  units = "in",
  dpi = 450
)

write_tsv(
  features |>
    select(analysis_id, isolate, study, orf, start, end, strand),
  file.path(output_dir, "densovirus_genome_organization_features.tsv")
)

write_tsv(
  sequences |>
    select(analysis_id, isolate, study, length_nt, unresolved_orfs),
  file.path(output_dir, "densovirus_genome_organization_sequences.tsv")
)

message("Densovirus genome-organization figure completed.")
message("Results: ", normalizePath(output_dir))
