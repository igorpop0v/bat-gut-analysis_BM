# Custom MaAsLin3 summary figures for all taxonomic ranks.

library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggnewscale)

source("workflows/16s/r/local_paths.R")

analysis_dir <- file.path(
  results_dir,
  "maaslin3",
  "all_ranks_primary"
)

results_file <- file.path(
  analysis_dir,
  "all_ranks_all_results.tsv"
)

taxonomy_file <- file.path(
  results_dir,
  "taxonomy",
  "asv_taxonomy_parsed_to_genus.tsv"
)

output_dir <- file.path(
  results_dir,
  "maaslin3",
  "custom_figures"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ranks <- c(
  "phylum",
  "class",
  "order",
  "family",
  "genus"
)

comparison_levels <- c(
  "ABH vs. H",
  "AAH vs. H"
)

results <- read_tsv(
  results_file,
  show_col_types = FALSE
)

taxonomy <- read_tsv(
  taxonomy_file,
  show_col_types = FALSE
)

# Family names used before unresolved GTDB genus identifiers.
placeholder_families <- taxonomy |>
  filter(
    genus %in% c(
      "BHER01",
      "JALENY01",
      "GCF-002002485"
    )
  ) |>
  distinct(genus, family) |>
  deframe()

# q >= 0.05 is white.
# Significant q values become progressively darker.
q_min <- 0.05
q_floor <- 0.000001

q_limits <- -log10(
  c(q_min, q_floor)
)

q_breaks <- -log10(
  c(
    0.05,
    0.01,
    0.001,
    0.0001,
    0.000001
  )
)

q_labels <- c(
  "≥0.05",
  "0.01",
  "0.001",
  "0.0001",
  "≤0.000001"
)

q_colour_positions <- scales::rescale(
  c(
    q_limits[1],
    2,
    4,
    q_limits[2]
  )
)

abundance_colours <- c(
  "#FFFFFF",
  "#E7CEE8",
  "#C06FBE",
  "#8E138A"
)

prevalence_colours <- c(
  "#FFFFFF",
  "#D7ECEA",
  "#69B8B2",
  "#168C87"
)

q_colour_value <- function(q) {
  ifelse(
    !is.na(q) & q < q_min,
    -log10(pmax(q, q_floor)),
    q_limits[1]
  )
}

# Use a true minus sign in numeric labels.
minus_labels <- function(x) {
  labels <- scales::label_number(
    accuracy = 0.1,
    trim = TRUE
  )(x)
  
  sub(
    "^-",
    "−",
    labels
  )
}

# Format taxon labels as in the composition heatmaps.
taxon_labels <- function(labels, rank_name) {
  lapply(labels, function(label) {
    clean_label <- gsub(
      "_",
      " ",
      label,
      fixed = TRUE
    )
    
    # Unresolved GTDB genus identifiers.
    if (
      rank_name == "genus" &&
      label %in% names(placeholder_families)
    ) {
      family_name <- gsub(
        "_",
        " ",
        placeholder_families[[label]],
        fixed = TRUE
      )
      
      return(
        bquote(
          italic(.(family_name)) ~ .(clean_label)
        )
      )
    }
    
    # Keep "Unclassified" in regular font.
    if (
      rank_name %in% c("family", "genus") &&
      startsWith(clean_label, "Unclassified ")
    ) {
      parent_name <- sub(
        "^Unclassified ",
        "",
        clean_label
      )
      
      return(
        bquote(
          "Unclassified " * italic(.(parent_name))
        )
      )
    }
    
    # Family and genus names are italicized.
    if (rank_name %in% c("family", "genus")) {
      return(
        bquote(
          italic(.(clean_label))
        )
      )
    }
    
    clean_label
  })
}

make_rank_figure <- function(rank_name) {
  
  rank_results <- results |>
    filter(
      taxonomic_rank == rank_name
    )
  
  # Show a taxon only when abundance or prevalence
  # is significant in at least one comparison.
  selected_taxa <- rank_results |>
    filter(
      metadata == "group",
      !is.na(qval_individual),
      qval_individual < 0.05
    ) |>
    pull(feature) |>
    unique()
  
  if (length(selected_taxa) == 0) {
    message(
      "No significant taxa at rank: ",
      rank_name
    )
    
    return(NULL)
  }
  
  # Order taxa by their smallest individual q-value.
  taxon_order <- rank_results |>
    filter(
      metadata == "group",
      feature %in% selected_taxa
    ) |>
    group_by(feature) |>
    summarise(
      best_individual_q = if (
        all(is.na(qval_individual))
      ) {
        Inf
      } else {
        min(
          qval_individual,
          na.rm = TRUE
        )
      },
      .groups = "drop"
    ) |>
    arrange(best_individual_q) |>
    pull(feature)
  
  taxon_order <- unique(
    c(
      taxon_order,
      setdiff(
        selected_taxa,
        taxon_order
      )
    )
  )
  
  group_data <- rank_results |>
    filter(
      metadata == "group",
      feature %in% selected_taxa,
      value %in% c(
        "active",
        "post_hibernation"
      )
    ) |>
    mutate(
      comparison = recode(
        value,
        active = "ABH vs. H",
        post_hibernation = "AAH vs. H"
      ),
      
      comparison = factor(
        comparison,
        levels = comparison_levels
      ),
      
      feature = factor(
        feature,
        levels = rev(taxon_order)
      ),
      
      # Retain the original MaAsLin3 coefficient direction:
      # ABH - H and AAH - H.
      plot_coefficient = coef,
      plot_null = null_hypothesis,
      
      confidence_low =
        plot_coefficient - 1.96 * stderr,
      
      confidence_high =
        plot_coefficient + 1.96 * stderr,
      
      q_colour =
        q_colour_value(qval_individual)
    )
  
  valid_data <- group_data |>
    filter(
      is.finite(plot_coefficient),
      is.finite(stderr)
    )
  
  null_lines <- group_data |>
    filter(
      is.finite(plot_null)
    ) |>
    distinct(
      comparison,
      model,
      plot_null
    )
  
  x_values <- c(
    valid_data$confidence_low,
    valid_data$confidence_high,
    null_lines$plot_null
  )
  
  x_limit <- ceiling(
    max(
      abs(x_values),
      na.rm = TRUE
    ) * 1.08
  )
  
  if (!is.finite(x_limit) || x_limit < 1) {
    x_limit <- 1
  }
  
  coefficient_plot <- ggplot() +
    
    # Null hypotheses for abundance and prevalence.
    geom_vline(
      data = null_lines,
      aes(
        xintercept = plot_null,
        linetype = model,
        colour = model
      ),
      linewidth = 0.55
    ) +
    
    # 95% confidence intervals.
    geom_errorbar(
      data = valid_data,
      aes(
        y = feature,
        xmin = confidence_low,
        xmax = confidence_high
      ),
      orientation = "y",
      width = 0.16,
      linewidth = 0.55,
      colour = "black"
    ) +
    
    # Abundance associations.
    geom_point(
      data = filter(
        valid_data,
        model == "abundance"
      ),
      aes(
        x = plot_coefficient,
        y = feature,
        fill = q_colour,
        shape = model
      ),
      colour = "black",
      stroke = 0.65,
      size = 3.2
    ) +
    
    scale_fill_gradientn(
      name = "Abundance q",
      colours = abundance_colours,
      values = q_colour_positions,
      limits = q_limits,
      breaks = q_breaks,
      labels = q_labels,
      oob = scales::squish,
      na.value = "white",
      guide = guide_colorbar(
        order = 4,
        barheight = grid::unit(36, "mm"),
        barwidth = grid::unit(5, "mm")
      )
    ) +
    
    ggnewscale::new_scale_fill() +
    
    # Prevalence associations.
    geom_point(
      data = filter(
        valid_data,
        model == "prevalence"
      ),
      aes(
        x = plot_coefficient,
        y = feature,
        fill = q_colour,
        shape = model
      ),
      colour = "black",
      stroke = 0.65,
      size = 3.4
    ) +
    
    scale_fill_gradientn(
      name = "Prevalence q",
      colours = prevalence_colours,
      values = q_colour_positions,
      limits = q_limits,
      breaks = q_breaks,
      labels = q_labels,
      oob = scales::squish,
      na.value = "white",
      guide = guide_colorbar(
        order = 3,
        barheight = grid::unit(36, "mm"),
        barwidth = grid::unit(5, "mm")
      )
    ) +
    
    scale_shape_manual(
      name = "Association",
      values = c(
        abundance = 21,
        prevalence = 24
      ),
      labels = c(
        abundance = "Abundance",
        prevalence = "Prevalence"
      ),
      guide = guide_legend(
        order = 2,
        override.aes = list(
          fill = "white",
          size = 3.2
        )
      )
    ) +
    
    scale_linetype_manual(
      name = "Null hypothesis",
      values = c(
        abundance = "dashed",
        prevalence = "solid"
      ),
      labels = c(
        abundance = "Abundance",
        prevalence = "Prevalence"
      ),
      guide = guide_legend(order = 1)
    ) +
    
    scale_colour_manual(
      name = "Null hypothesis",
      values = c(
        abundance = "grey55",
        prevalence = "black"
      ),
      labels = c(
        abundance = "Abundance",
        prevalence = "Prevalence"
      ),
      guide = guide_legend(order = 1)
    ) +
    
    scale_x_continuous(
      limits = c(
        -x_limit,
        x_limit
      ),
      labels = minus_labels,
      expand = expansion(mult = 0.03)
    ) +
    
    scale_y_discrete(
      labels = function(x) {
        taxon_labels(
          x,
          rank_name
        )
      }
    ) +
    
    facet_grid(
      . ~ comparison,
      drop = FALSE
    ) +
    
    labs(
      x = "β coefficient",
      y = "Taxon"
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      panel.grid.minor = element_blank(),
      
      panel.grid.major.y = element_line(
        colour = "grey90",
        linewidth = 0.35
      ),
      
      panel.spacing.x = grid::unit(
        3,
        "mm"
      ),
      
      strip.background = element_rect(
        fill = "white",
        colour = "black"
      ),
      
      strip.text = element_text(
        size = 13,
        face = "bold"
      ),
      
      axis.text.x = element_text(
        size = 11,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        size = 11.5,
        colour = "black"
      ),
      
      axis.title = element_text(
        size = 14
      ),
      
      legend.position = "right",
      
      legend.title = element_text(
        size = 14,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      ),
      
      legend.box = "vertical",
      
      plot.margin = margin(
        8,
        8,
        8,
        8
      )
    )
  
  figure_height <- max(
    7,
    2.8 + 0.30 * length(taxon_order)
  )
  
  figure_width <- 11
  
  # Save the exact data used in each figure.
  write_tsv(
    group_data |>
      mutate(
        feature = as.character(feature)
      ),
    file.path(
      output_dir,
      paste0(
        rank_name,
        "_plot_data.tsv"
      )
    )
  )
  
  # PNG for quick viewing.
  ggsave(
    file.path(
      output_dir,
      paste0(
        "maaslin3_",
        rank_name,
        "_custom.png"
      )
    ),
    coefficient_plot,
    width = figure_width,
    height = figure_height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  # Vector PDF for the thesis or publication.
  ggsave(
    file.path(
      output_dir,
      paste0(
        "maaslin3_",
        rank_name,
        "_custom.pdf"
      )
    ),
    coefficient_plot,
    device = cairo_pdf,
    width = figure_width,
    height = figure_height,
    units = "in",
    bg = "white"
  )
  
  tibble(
    taxonomic_rank = rank_name,
    plotted_taxa = length(taxon_order),
    figure_width_inches = figure_width,
    figure_height_inches = figure_height
  )
}

figure_summary <- lapply(
  ranks,
  make_rank_figure
) |>
  bind_rows()

write_tsv(
  figure_summary,
  file.path(
    output_dir,
    "custom_figure_summary.tsv"
  )
)

message(
  "Custom MaAsLin3 figures completed."
)

message(
  "Output directory: ",
  output_dir
)