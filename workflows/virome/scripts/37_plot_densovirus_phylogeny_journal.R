#!/usr/bin/env Rscript

# Plot the study-inclusive NSP1 + SP1 densovirus phylogeny in the visual
# style used in the previous BatShotMetaFlow publication.

required_packages <- c(
  "ape",
  "dplyr",
  "ggimage",
  "ggnewscale",
  "ggplot2",
  "ggtree",
  "rsvg",
  "viridis"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the missing R packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggimage)
  library(ggnewscale)
  library(ggplot2)
  library(ggtree)
  library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    "Usage: 37_plot_densovirus_phylogeny_journal.R ",
    "TREEFILE TREE_METADATA.tsv OUTPUT_DIR FLAGS_DIR"
  )
}

tree_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]
flags_dir <- args[4]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(flags_dir, recursive = TRUE, showWarnings = FALSE)

minimum_alrt <- 80
minimum_ufboot <- 95

outgroup_ids <- c(
  "NC_011545.2",
  "NC_077018.1",
  "NC_015115.1",
  "NC_011317.1",
  "NC_012636.1"
)

group_colors <- c(
  "ABH" = "#0000FF",
  "H" = "#FF0000",
  "AAH" = "#BDBDBD"
)

study_colors <- c(
  "Popov et al., 2025" = "#E69F00",
  "Yang et al., 2016" = "#333333",
  "Armién et al., 2023" = "#009E73",
  "This study: ABH" = group_colors[["ABH"]],
  "This study: H" = group_colors[["H"]],
  "This study: AAH" = group_colors[["AAH"]]
)

metadata <- read.delim(
  metadata_file,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Name",
  "Full.Name",
  "Country",
  "Year",
  "Host",
  "Host.icon.query",
  "flag_code",
  "Dataset",
  "Group",
  "Study"
)

if (!all(required_columns %in% colnames(metadata))) {
  stop("Tree metadata does not contain all required columns.")
}

tree <- read.tree(tree_file)

if (!all(outgroup_ids %in% tree$tip.label)) {
  stop("One or more outgroup sequences are absent from the tree.")
}

missing_metadata <- setdiff(tree$tip.label, metadata$Name)

if (length(missing_metadata) > 0) {
  stop(
    "No metadata for these tree tips: ",
    paste(missing_metadata, collapse = ", ")
  )
}

tree <- root(
  tree,
  outgroup = outgroup_ids,
  resolve.root = TRUE
)

tree <- ladderize(tree, right = FALSE)

metadata <- metadata %>%
  filter(Name %in% tree$tip.label) %>%
  mutate(
    Year = na_if(Year, "ND"),
    Year = factor(
      Year,
      levels = sort(unique(Year[!is.na(Year)]))
    ),
    Group = na_if(Group, "ND"),
    Group = factor(Group, levels = c("ABH", "H", "AAH")),
    Study = na_if(Study, "ND"),
    Study = factor(Study, levels = names(study_colors)),
    flag_code = na_if(flag_code, "ND"),
    Host.icon.query = na_if(Host.icon.query, "ND")
  )

# Download circular country flags once and reuse the local copies.
flag_codes <- sort(unique(metadata$flag_code[!is.na(metadata$flag_code)]))

for (flag_code in flag_codes) {
  flag_file <- file.path(flags_dir, paste0(flag_code, ".svg"))

  if (!file.exists(flag_file)) {
    flag_url <- paste0(
      "https://raw.githubusercontent.com/",
      "HatScripts/circle-flags/gh-pages/flags/",
      flag_code,
      ".svg"
    )

    download.file(
      flag_url,
      flag_file,
      mode = "wb",
      quiet = TRUE
    )
  }
}

metadata$flag_path <- ifelse(
  is.na(metadata$flag_code),
  NA_character_,
  file.path(flags_dir, paste0(metadata$flag_code, ".svg"))
)

# Query PhyloPic once. The resolved identifiers are saved for reproducibility.
phylopic_file <- file.path(
  output_dir,
  "densovirus_host_phylopic_uids.tsv"
)

host_queries <- sort(unique(
  metadata$Host.icon.query[!is.na(metadata$Host.icon.query)]
))

if (file.exists(phylopic_file)) {
  phylopic_uids <- read.delim(
    phylopic_file,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  )
} else {
  resolve_phylopic_uid <- function(host_query) {
    tryCatch(
      ggimage::phylopic_uid(
        host_query,
        seed = 123
      )$uid[[1]],
      error = function(error) {
        warning(
          "No PhyloPic image for '",
          host_query,
          "'; this host will be plotted without a silhouette."
        )
        NA_character_
      }
    )
  }

  phylopic_uids <- data.frame(
    name = host_queries,
    uid = vapply(
      host_queries,
      resolve_phylopic_uid,
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  write.table(
    phylopic_uids,
    phylopic_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

metadata <- metadata %>%
  left_join(
    phylopic_uids,
    by = c("Host.icon.query" = "name")
  )

if (any(is.na(metadata$uid))) {
  warning(
    "PhyloPic was not resolved for: ",
    paste(
      unique(metadata$Host[is.na(metadata$uid)]),
      collapse = ", "
    )
  )
}

# Cache the resolved PhyloPic silhouettes locally. This keeps the plotting
# layers independent from the PhyloPic API and makes later reruns offline-safe.
host_icons_dir <- file.path(output_dir, "host_icons")
dir.create(host_icons_dir, recursive = TRUE, showWarnings = FALSE)

download_host_icon <- function(uid) {
  if (is.na(uid) || uid == "") {
    return(NA_character_)
  }

  host_icon_file <- file.path(host_icons_dir, paste0(uid, ".svg"))

  if (!file.exists(host_icon_file)) {
    host_icon_url <- paste0(
      "https://images.phylopic.org/images/",
      uid,
      "/vector.svg"
    )

    tryCatch(
      download.file(
        host_icon_url,
        host_icon_file,
        mode = "wb",
        quiet = TRUE
      ),
      error = function(error) {
        warning(
          "Could not download PhyloPic image ",
          uid,
          ": ",
          conditionMessage(error)
        )
      }
    )
  }

  if (file.exists(host_icon_file)) {
    host_icon_file
  } else {
    NA_character_
  }
}

metadata$host_icon_path <- vapply(
  metadata$uid,
  download_host_icon,
  character(1)
)

write.table(
  metadata,
  file.path(output_dir, "densovirus_tree_metadata_resolved.tsv"),
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

      parts <- strsplit(value, "/", fixed = TRUE)[[1]]

      if (length(parts) < position) {
        return(NA_real_)
      }

      suppressWarnings(as.numeric(parts[position]))
    },
    numeric(1)
  )
}

# Choose a conventional, easy-to-read scale-bar length close to one fifth of
# the observed tree width (for example, 1 or 0.01 substitutions per site).
choose_scale_bar_width <- function(tree_width) {
  target_width <- tree_width / 5
  magnitude <- 10^floor(log10(target_width))
  scaled_target <- target_width / magnitude

  nice_multiplier <- if (scaled_target >= 5) {
    5
  } else if (scaled_target >= 2) {
    2
  } else {
    1
  }

  nice_multiplier * magnitude
}

make_tree_plot <- function(
    phy,
    plot_metadata,
    plot_type = c("full", "focal")) {

  plot_type <- match.arg(plot_type)
  phy <- ladderize(phy, right = FALSE)

  attached_metadata <- plot_metadata %>%
    filter(Name %in% phy$tip.label) %>%
    transmute(
      label = Name,
      Full.Name,
      Country,
      Year,
      Host,
      host_icon_path,
      flag_path,
      Group,
      Study
    )

  meta_year <- attached_metadata %>%
    select(label, Year) %>%
    tibble::column_to_rownames("label")

  meta_study <- attached_metadata %>%
    select(label, Study) %>%
    tibble::column_to_rownames("label")

  base_plot <- ggtree(
    phy,
    linewidth = ifelse(plot_type == "full", 0.32, 0.48)
  )

  base_plot$data$sh_alrt <- parse_support(base_plot$data$label, 1)
  base_plot$data$ufboot <- parse_support(base_plot$data$label, 2)

  branch_data <- base_plot$data %>%
    filter(!isTip) %>%
    mutate(
      strong_support = case_when(
        is.na(sh_alrt) & is.na(ufboot) ~ TRUE,
        !is.na(sh_alrt) & !is.na(ufboot) &
          sh_alrt >= minimum_alrt &
          ufboot >= minimum_ufboot ~ TRUE,
        TRUE ~ FALSE
      )
    )

  maximum_x <- max(base_plot$data$x, na.rm = TRUE)
  scale_bar_width <- choose_scale_bar_width(maximum_x)

  if (plot_type == "full") {
    label_size <- 2.35
    host_size <- 0.021
    flag_size <- 0.019
    host_offset <- maximum_x * 0.47
    country_offset <- maximum_x * 0.57
    year_offset <- maximum_x * 0.67
    x_limit <- maximum_x * 1.82
    header_size <- 3.2
  } else {
    label_size <- 3.2
    host_size <- 0.037
    flag_size <- 0.035
    host_offset <- maximum_x * 1.35
    country_offset <- maximum_x * 1.55
    year_offset <- maximum_x * 1.75
    study_offset <- maximum_x * 1.93
    x_limit <- maximum_x * 3.15
    header_size <- 3.8
  }

  tip_annotations <- base_plot$data %>%
    filter(isTip) %>%
    select(label, x, y) %>%
    left_join(attached_metadata, by = "label") %>%
    mutate(
      host_x = maximum_x + host_offset,
      country_x = maximum_x + country_offset
    )

  host_annotations <- tip_annotations %>%
    filter(!is.na(host_icon_path), host_icon_path != "")

  flag_annotations <- tip_annotations %>%
    filter(!is.na(flag_path), flag_path != "")

  p <- base_plot %<+% attached_metadata

  p <- p +
    ggnewscale::new_scale_color() +
    geom_tree(
      data = branch_data,
      aes(color = strong_support),
      linewidth = ifelse(plot_type == "full", 0.32, 0.48)
    ) +
    scale_color_manual(
      values = c("TRUE" = "black", "FALSE" = "grey75"),
      guide = "none"
    ) +
    geom_segment(
      data = tip_annotations,
      aes(
        x = maximum_x,
        xend = host_x,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      color = "grey60",
      linetype = "dotted",
      linewidth = 0.25
    ) +
    geom_label(
      data = tip_annotations,
      aes(x = x, y = y, label = Full.Name),
      inherit.aes = FALSE,
      hjust = 0,
      color = "black",
      fill = "white",
      linewidth = 0,
      label.padding = grid::unit(0.02, "lines"),
      size = label_size
    ) +
    ggimage::geom_image(
      data = host_annotations,
      aes(x = host_x, y = y, image = host_icon_path),
      inherit.aes = FALSE,
      size = host_size,
      na.rm = TRUE
    ) +
    ggimage::geom_image(
      data = flag_annotations,
      aes(x = country_x, y = y, image = flag_path),
      inherit.aes = FALSE,
      size = flag_size,
      na.rm = TRUE
    )

  p <- suppressMessages({
    year_plot <- gheatmap(
      p,
      meta_year,
      width = 0.045,
      offset = year_offset,
      color = "black",
      font.size = header_size,
      colnames_offset_y = -0.05
    )

    year_plot +
      scale_fill_viridis_d(
        option = "D",
        name = "Year",
        na.value = "white",
        na.translate = TRUE
      )
  })

  if (plot_type == "focal") {
    p <- p + ggnewscale::new_scale_fill()

    p <- suppressMessages({
      study_plot <- gheatmap(
        p,
        meta_study,
        width = 0.045,
        offset = study_offset,
        color = "black",
        font.size = header_size,
        colnames_offset_y = -0.05
      )

      study_plot +
        scale_fill_manual(
          values = study_colors,
          breaks = names(study_colors),
          name = "Study",
          na.value = "white",
          na.translate = FALSE,
          drop = FALSE
        )
    })
  }

  p +
    geom_treescale(
      x = 0,
      y = -0.65,
      width = scale_bar_width,
      offset = 0.22,
      linesize = ifelse(plot_type == "full", 0.45, 0.6),
      fontsize = ifelse(plot_type == "full", 2.8, 3.4),
      family = "sans"
    ) +
    annotate(
      "text",
      x = maximum_x + host_offset,
      y = 0,
      label = "Host",
      size = header_size
    ) +
    annotate(
      "text",
      x = maximum_x + country_offset,
      y = 0,
      label = "Country",
      size = header_size
    ) +
    xlim(0, x_limit) +
    theme_tree() +
    theme(
      text = element_text(
        family = "sans",
        color = "black",
        size = ifelse(plot_type == "full", 9, 11)
      ),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "black"),
      plot.margin = margin(8, 16, 12, 8)
    )
}

# Complete tree.
full_tree_plot <- make_tree_plot(
  tree,
  metadata,
  plot_type = "full"
)

# Focused view containing all study sequences and their two closest references.
study_ids <- grep("^(Current|Previous)_", tree$tip.label, value = TRUE)

focal_ids <- intersect(
  c(study_ids, "NC_031450.1", "MW628494.1"),
  tree$tip.label
)

focal_tree <- keep.tip(tree, focal_ids)

focal_tree_plot <- make_tree_plot(
  focal_tree,
  metadata,
  plot_type = "focal"
)

ggsave(
  file.path(output_dir, "densovirus_partitioned_full_tree_journal.png"),
  full_tree_plot,
  width = 15,
  height = 18,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(output_dir, "densovirus_partitioned_full_tree_journal.pdf"),
  full_tree_plot,
  width = 15,
  height = 18,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  file.path(output_dir, "densovirus_partitioned_focal_tree_journal.png"),
  focal_tree_plot,
  width = 14,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(output_dir, "densovirus_partitioned_focal_tree_journal.pdf"),
  focal_tree_plot,
  width = 14,
  height = 9,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

cat("Journal-style densovirus trees created.\n")
cat("Results:", output_dir, "\n")
