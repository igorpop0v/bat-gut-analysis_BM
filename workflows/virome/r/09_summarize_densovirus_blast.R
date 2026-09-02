# Summarize BLASTn HSPs for densovirus candidate-reference pairs.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage: Rscript 09_summarize_densovirus_blast.R",
      "<all_hits.tsv> <output_directory> <candidate_fasta>"
    )
  )
}

all_hits_file <- args[1]
output_dir <- args[2]
candidate_fasta <- args[3]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setlocale("LC_NUMERIC", "C")
options(OutDec = ".")

column_names <- c(
  "query_id",
  "reference_id",
  "identity_pct",
  "alignment_length",
  "mismatches",
  "gap_openings",
  "query_start",
  "query_end",
  "query_length",
  "reference_start",
  "reference_end",
  "reference_length",
  "e_value",
  "bitscore"
)

hits <- read.delim(
  all_hits_file,
  header = FALSE,
  sep = "\t",
  col.names = column_names,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(hits) == 0) {
  stop("No BLASTn hits were found.")
}

# Calculate the total number of uniquely covered positions.
merge_interval_length <- function(starts, ends) {

  interval_table <- data.frame(
    start = pmin(starts, ends),
    end = pmax(starts, ends)
  )

  interval_table <- interval_table[
    order(interval_table$start, interval_table$end),
    ,
    drop = FALSE
  ]

  current_start <- interval_table$start[1]
  current_end <- interval_table$end[1]
  covered_length <- 0

  if (nrow(interval_table) > 1) {

    for (i in 2:nrow(interval_table)) {

      if (interval_table$start[i] <= current_end + 1) {

        current_end <- max(
          current_end,
          interval_table$end[i]
        )

      } else {

        covered_length <- covered_length +
          current_end - current_start + 1

        current_start <- interval_table$start[i]
        current_end <- interval_table$end[i]
      }
    }
  }

  covered_length + current_end - current_start + 1
}

pair_key <- paste(
  hits$query_id,
  hits$reference_id,
  sep = "\034"
)

split_hits <- split(
  hits,
  pair_key
)

pair_summary <- lapply(
  split_hits,
  function(current_hits) {

    query_covered_bases <- merge_interval_length(
      current_hits$query_start,
      current_hits$query_end
    )

    reference_covered_bases <- merge_interval_length(
      current_hits$reference_start,
      current_hits$reference_end
    )

    best_hsp <- current_hits[
      which.max(current_hits$bitscore),
      ,
      drop = FALSE
    ]

    data.frame(
      query_id = current_hits$query_id[1],
      reference_id = current_hits$reference_id[1],
      best_hsp_identity_pct = best_hsp$identity_pct[1],
      weighted_hsp_identity_pct = weighted.mean(
        current_hits$identity_pct,
        current_hits$alignment_length
      ),
      hsp_count = nrow(current_hits),
      query_length = max(current_hits$query_length),
      reference_length = max(current_hits$reference_length),
      query_covered_bases = query_covered_bases,
      reference_covered_bases = reference_covered_bases,
      query_coverage_pct = 100 *
        query_covered_bases /
        max(current_hits$query_length),
      reference_coverage_pct = 100 *
        reference_covered_bases /
        max(current_hits$reference_length),
      minimum_e_value = min(current_hits$e_value),
      maximum_bitscore = max(current_hits$bitscore),
      total_bitscore = sum(current_hits$bitscore),
      stringsAsFactors = FALSE
    )
  }
)

pair_summary <- do.call(
  rbind,
  pair_summary
)

row.names(pair_summary) <- NULL

# Select the best reference using reference coverage first.
pair_summary <- pair_summary[
  order(
    pair_summary$query_id,
    -pair_summary$reference_coverage_pct,
    -pair_summary$query_coverage_pct,
    -pair_summary$maximum_bitscore
  ),
  ,
  drop = FALSE
]

best_hits <- pair_summary[
  !duplicated(pair_summary$query_id),
  ,
  drop = FALSE
]

round_columns <- c(
  "best_hsp_identity_pct",
  "weighted_hsp_identity_pct",
  "query_coverage_pct",
  "reference_coverage_pct"
)

pair_summary[round_columns] <- lapply(
  pair_summary[round_columns],
  round,
  digits = 3
)

best_hits[round_columns] <- lapply(
  best_hits[round_columns],
  round,
  digits = 3
)

write.table(
  pair_summary,
  file.path(
    output_dir,
    "densovirus_blastn_pair_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  best_hits,
  file.path(
    output_dir,
    "densovirus_blastn_best_hits.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Identify candidate contigs without matches to previous genomes.
fasta_headers <- readLines(candidate_fasta)
fasta_headers <- fasta_headers[startsWith(fasta_headers, ">")]
candidate_ids <- sub("^>", "", fasta_headers)
candidate_ids <- sub("[[:space:]].*$", "", candidate_ids)

candidates_without_hits <- setdiff(
  candidate_ids,
  unique(hits$query_id)
)

write.table(
  data.frame(query_id = candidates_without_hits),
  file.path(
    output_dir,
    "densovirus_candidates_without_previous_reference_hit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("BLASTn HSP summarization completed.")
message("Candidate-reference pairs: ", nrow(pair_summary))
message("Candidates with previous-genome hits: ", nrow(best_hits))
message("Candidates without previous-genome hits: ",
        length(candidates_without_hits))
