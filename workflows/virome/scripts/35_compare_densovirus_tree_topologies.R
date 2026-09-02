#!/usr/bin/env Rscript

library(ape)

args <- commandArgs(trailingOnly = TRUE)

tree_dir <- args[1]
output_dir <- args[2]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

alrt_threshold <- 80
ufboot_threshold <- 95

tree_files <- c(
  primary = file.path(
    tree_dir,
    "primary",
    "densovirus_primary_nsp1_sp1.treefile"
  ),
  study_inclusive = file.path(
    tree_dir,
    "study_inclusive",
    "densovirus_study_inclusive_nsp1_sp1.treefile"
  ),
  sensitivity = file.path(
    tree_dir,
    "sensitivity",
    "densovirus_sensitivity_nsp1_sp1.treefile"
  )
)

if (!all(file.exists(tree_files))) {
  stop("One or more IQ-TREE tree files are missing.")
}

unroot_tree <- function(tree) {

  if (is.rooted(tree)) {
    tree <- unroot(tree)
  }

  tree
}

collapse_weak_nodes <- function(
  tree,
  minimum_alrt = 80,
  minimum_ufboot = 95
) {

  number_of_tips <- Ntip(tree)

  for (index in seq_along(tree$node.label)) {

    node <- number_of_tips + index

    support <- strsplit(
      tree$node.label[index],
      "/",
      fixed = TRUE
    )[[1]]

    alrt <- suppressWarnings(as.numeric(support[1]))
    ufboot <- suppressWarnings(as.numeric(support[2]))

    strong_node <- (
      length(support) >= 2 &&
      !is.na(alrt) &&
      !is.na(ufboot) &&
      alrt >= minimum_alrt &&
      ufboot >= minimum_ufboot
    )

    parent_edge <- which(tree$edge[, 2] == node)

    if (!strong_node && length(parent_edge) == 1) {
      tree$edge.length[parent_edge] <- 0
    }
  }

  di2multi(tree, tol = 1e-10)
}

calculate_rf <- function(tree_1, tree_2, tips) {

  tree_1 <- keep.tip(tree_1, tips)
  tree_2 <- keep.tip(tree_2, tips)

  tree_1 <- unroot_tree(tree_1)
  tree_2 <- unroot_tree(tree_2)

  rf_distance <- as.numeric(
    dist.topo(
      tree_1,
      tree_2,
      method = "PH85"
    )
  )

  maximum_rf <- 2 * (length(tips) - 3)

  data.frame(
    tips = length(tips),
    rf_distance = rf_distance,
    normalized_rf = rf_distance / maximum_rf
  )
}

raw_trees <- lapply(
  tree_files,
  read.tree
)

robust_trees <- lapply(
  raw_trees,
  collapse_weak_nodes,
  minimum_alrt = alrt_threshold,
  minimum_ufboot = ufboot_threshold
)

tree_pairs <- combn(
  names(raw_trees),
  2,
  simplify = FALSE
)

results <- lapply(
  tree_pairs,
  function(pair) {

    common_tips <- intersect(
      raw_trees[[pair[1]]]$tip.label,
      raw_trees[[pair[2]]]$tip.label
    )

    study_tips <- grep(
      "^(Current|Previous)_",
      common_tips,
      value = TRUE
    )

    raw_all <- calculate_rf(
      raw_trees[[pair[1]]],
      raw_trees[[pair[2]]],
      common_tips
    )

    robust_all <- calculate_rf(
      robust_trees[[pair[1]]],
      robust_trees[[pair[2]]],
      common_tips
    )

    raw_study <- calculate_rf(
      raw_trees[[pair[1]]],
      raw_trees[[pair[2]]],
      study_tips
    )

    robust_study <- calculate_rf(
      robust_trees[[pair[1]]],
      robust_trees[[pair[2]]],
      study_tips
    )

    data.frame(
      tree_1 = pair[1],
      tree_2 = pair[2],
      common_tips = raw_all$tips,
      raw_rf = raw_all$rf_distance,
      raw_normalized_rf = raw_all$normalized_rf,
      robust_rf = robust_all$rf_distance,
      robust_normalized_rf = robust_all$normalized_rf,
      study_tips = raw_study$tips,
      study_raw_rf = raw_study$rf_distance,
      study_raw_normalized_rf = raw_study$normalized_rf,
      study_robust_rf = robust_study$rf_distance,
      study_robust_normalized_rf =
        robust_study$normalized_rf
    )
  }
)

results <- do.call(
  rbind,
  results
)

write.table(
  results,
  file.path(
    output_dir,
    "densovirus_tree_topology_comparison.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Tree topology comparison completed.\n")
cat("Results:", output_dir, "\n")

print(results)
