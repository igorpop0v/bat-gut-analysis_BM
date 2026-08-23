# Parse GTDB taxonomy and fill unresolved ranks up to genus level.

library(qiime2R)
library(dplyr)
library(tibble)
library(readr)

# Local input and output paths.
source("workflows/16s/r/local_paths.R")

taxonomy_file <- file.path(
  input_dir,
  "asv_taxonomy_gtdb_r232.qza"
)

asv_counts_file <- file.path(
  input_dir,
  "asv_counts_unrarefied.tsv"
)

taxonomy_output_dir <- file.path(
  results_dir,
  "taxonomy"
)

output_file <- file.path(
  taxonomy_output_dir,
  "asv_taxonomy_parsed_to_genus.tsv"
)

dir.create(
  taxonomy_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Read the QIIME 2 taxonomy artifact.
taxonomy_raw <- read_qza(taxonomy_file)$data

# Parse the taxonomy string into separate ranks.
taxonomy_parsed <- parse_taxonomy(taxonomy_raw) |>
  rownames_to_column("feature_id") |>
  rename(
    domain = Kingdom,
    phylum = Phylum,
    class = Class,
    order = Order,
    family = Family,
    genus = Genus
  ) |>
  select(
    feature_id,
    domain,
    phylum,
    class,
    order,
    family,
    genus
  )

# Preserve the confidence value from the QIIME 2 classifier.
taxonomy_confidence <- taxonomy_raw |>
  transmute(
    feature_id = Feature.ID,
    confidence = Confidence
  )

taxonomy_clean <- taxonomy_parsed |>
  left_join(
    taxonomy_confidence,
    by = "feature_id"
  )

# Remove remaining rank prefixes, including the GTDB domain prefix d__.
rank_columns <- c(
  "domain",
  "phylum",
  "class",
  "order",
  "family",
  "genus"
)

taxonomy_clean[rank_columns] <- lapply(
  taxonomy_clean[rank_columns],
  function(x) {
    x <- trimws(as.character(x))
    x <- sub("^[dkpcofgs]__", "", x)
    x[x == ""] <- NA_character_
    x
  }
)

# Identify missing or unresolved rank names.
is_unresolved <- function(x) {
  x_lower <- tolower(trimws(x))

  is.na(x) |
    x_lower == "" |
    grepl(
      "^(uncultured|unclassified|unknown|unidentified)",
      x_lower
    )
}

# Fill a missing domain if one is encountered.
taxonomy_clean$domain[
  is_unresolved(taxonomy_clean$domain)
] <- "Unclassified"

# Fill each missing rank using the nearest classified parent rank.
for (rank_index in 2:length(rank_columns)) {

  current_rank <- rank_columns[rank_index]
  parent_rank <- rank_columns[rank_index - 1]

  unresolved_rows <- is_unresolved(
    taxonomy_clean[[current_rank]]
  )

  parent_values <- taxonomy_clean[[parent_rank]][
    unresolved_rows
  ]

  replacement_values <- ifelse(
    startsWith(parent_values, "Unclassified_") |
      parent_values == "Unclassified",
    parent_values,
    paste0("Unclassified_", parent_values)
  )

  taxonomy_clean[[current_rank]][
    unresolved_rows
  ] <- replacement_values
}

# Read ASV IDs from the biological-sample count table.
asv_ids <- read_tsv(
  asv_counts_file,
  skip = 1,
  show_col_types = FALSE,
  name_repair = "minimal"
)[[1]]

# Retain only ASVs present in the biological-sample count table
# and arrange taxonomy in the same order.
taxonomy_clean <- tibble(
  feature_id = asv_ids
) |>
  left_join(
    taxonomy_clean,
    by = "feature_id"
  )

if (anyNA(taxonomy_clean$domain)) {
  stop("Some ASVs in the count table have no taxonomy assignment.")
}

# Write a tab-separated file without row names.
write_tsv(
  taxonomy_clean,
  output_file,
  na = ""
)

genus_unclassified <- startsWith(
  taxonomy_clean$genus,
  "Unclassified"
)

message("Taxonomy parsing finished.")
message("Biological ASVs: ", nrow(taxonomy_clean))
message(
  "ASVs classified to genus: ",
  sum(!genus_unclassified)
)
message(
  "ASVs unresolved at genus: ",
  sum(genus_unclassified)
)
message("Output: ", output_file)
