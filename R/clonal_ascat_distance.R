####################################################################################################
#' This function computes various "distances", which are used as penalties for a copy number solution
#' One such distance is an estimate of the proportion of the tumour genome which is clonal.
#' For each segment of the genome, we test the null hypothesis is that
#' the tumour genome segment in question is "clonal". The alternative hypothesis is that
#' the tumour genome segment in question exhibits "sub-clonal" variation.
#' Cleaned version of clonal distance calculation
#' @noRd
calc_distance_clonal <- function(
  segs, dist_choice, rho, psi, gamma_param, read_depth,
  siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR,
  uninformative_baf_threshold
) {
  # Filter informative segments up front
  s <- segs[segs[, "b"] > uninformative_baf_threshold, , drop = FALSE]

  # Handle empty case immediately to match original logic
  if (nrow(s) == 0) {
    return(list(
      distance_value = 0,
      minimise = FALSE,
      max_clonal_segment = 0,
      ref_maj = NA,
      ref_min = NA
    ))
  }

  # Map: Calculate metrics for every segment
  stats <- lapply(seq_len(nrow(s)), function(i) {
    row <- s[i, ]

    seg_info <- is_segment_clonal(
      row["r"], row["b"], row["length"], row["size"],
      row["mean"], row["sd"], rho, psi, gamma_param,
      siglevel_BAF, maxdist_BAF
    )

    err_info <- calc_standardised_error(
      row["r"], row["b"], row["length"], row["size"],
      row["mean"], row["sd"], rho, psi, gamma_param, maxdist_BAF
    )

    ln_lratio <- calc_ln_likelihood_ratio(
      row["r"], row["b"], row["length"], row["size"],
      row["mean"], read_depth, rho, psi, gamma_param, maxdist_BAF
    )

    list(
      is_clonal   = seg_info$is_clonal,
      is_balanced = seg_info$balanced,
      nMaj        = seg_info$nMaj,
      nMin        = seg_info$nMin,
      tvar_sq     = err_info$tvar^2,
      included    = err_info$included_segment,
      ln_lratio   = ln_lratio,
      size        = row["length"],
      b_diff_sq   = (row["b"] - row["mean"])^2
    )
  })

  # Aggregate Data
  sizes <- sapply(stats, `[[`, "size")
  is_clonal <- sapply(stats, `[[`, "is_clonal")
  is_balanced <- sapply(stats, `[[`, "is_balanced")

  genome_size <- sum(sizes)
  clonal_genome_size <- sum(sizes[is_clonal])
  n_inc <- sum(sapply(stats, `[[`, "included"))
  total_segs <- nrow(s)

  # Handle "Max Clonal Segment" logic (DCW 160314 balanced check)
  max_idx <- 0
  ref_maj <- NA
  ref_min <- NA
  potential_indices <- which(is_clonal & !is_balanced)

  if (length(potential_indices) > 0) {
    best_match_idx <- potential_indices[which.max(sizes[potential_indices])]
    max_idx <- best_match_idx
    ref_maj <- stats[[best_match_idx]]$nMaj
    ref_min <- stats[[best_match_idx]]$nMin
  }

  # Compute Final Distance with pmax safety checks
  res <- switch(as.character(dist_choice),
    "0" = list(v = clonal_genome_size / pmax(genome_size, 1e-10), m = FALSE), # Clonal Prop
    "1" = list(v = sum(sapply(stats, `[[`, "tvar_sq")) / pmax(n_inc, 1), m = TRUE),
    "2" = list(v = sum(sapply(stats, `[[`, "b_diff_sq")) / pmax(total_segs, 1), m = TRUE),
    "3" = list(v = sum(sapply(stats, function(x) x$size * x$b_diff_sq)) / pmax(genome_size, 1e-10), m = TRUE),
    "4" = list(v = sum(sapply(stats, `[[`, "ln_lratio")), m = FALSE)
  )

  return(list(
    distance_value     = res$v,
    minimise           = res$m,
    max_clonal_segment = max_idx,
    ref_maj            = ref_maj,
    ref_min            = ref_min
  ))
}

####################################################################################################
#' This function computes various "distances", which are used as penalties for a copy number solution.
#' This function is called when searching for a clonal copy number solution.
#' One such distance is an estimate of the proportion of the tumour genome which is clonal.
#' For each segment of the genome, we test the null hypothesis is that
#' the tumour genome segment in question is "clonal". The alternative hypothesis is that
#' the tumour genome segment in question exhibits "sub-clonal" variation.
#' @noRd
calc_distance <- function(
  segs, dist_choice, rho, psi, gamma_param,
  uninformative_baf_threshold = 0.51
) {
  s <- segs

  # Shared calculation for nA and nB (identical for all choices)
  # Pre-calculate common term to keep it clean
  common_multiplier <- 2^(s[, "r"] / gamma_param) * ((1 - rho) * 2 + rho * psi)
  nA <- (rho - 1 - (s[, "b"] - 1) * common_multiplier) / rho
  nB <- (rho - 1 + s[, "b"] * common_multiplier) / rho

  # Identify Minor and Major alleles
  # We compare sums once to determine the assignment
  if (sum(nA, na.rm = TRUE) < sum(nB, na.rm = TRUE)) {
    nMinor <- nA
    nMajor <- nB
  } else {
    nMinor <- nB
    nMajor <- nA
  }
  # Helper function for the "roundness" penalty used in all choices
  get_penalty <- function(n) (abs(n - pmax(round(n), 0)))

  # Specific distance logic
  if (dist_choice == 0) {
    # Original ASCAT distance
    weights <- ifelse(s[, "b"] <= uninformative_baf_threshold, 0.05, 1)
    dist_value <- sum(get_penalty(nMinor)^2 * s[, "length"] * weights, na.rm = TRUE)
    minimise <- TRUE
  } else {
    # All choices 1, 2, and 3 use (0.5 - penalty)^2
    pMinor <- (0.5 - get_penalty(nMinor))^2
    pMajor <- (0.5 - get_penalty(nMajor))^2
    minimise <- FALSE

    if (dist_choice == 1) {
      dist_value <- sum(pMinor * s[, "length"], na.rm = TRUE)
    } else if (dist_choice == 2) {
      dist_value <- 0.5 * sum((pMinor + pMajor) * s[, "length"], na.rm = TRUE)
    } else if (dist_choice == 3) {
      # Penalty for homozygous deletions
      hom_del <- nMinor < 0.5 & nMajor < 0.5 & nMinor >= 0 & nMajor >= 0

      # Multiply the penalty by 4 and the length by 2 for hom_dels
      segs_penalty <- (pMinor + pMajor)
      segs_penalty[hom_del] <- segs_penalty[hom_del] * 4
      dist_value <- 0.5 * sum(segs_penalty * (s[, "length"] * ifelse(hom_del, 2, 1)), na.rm = TRUE)
    }
  }
  return(list(distance_value = dist_value, minimise = minimise))
}


####################################################################################################
#' function to create the distance matrix (distance for a range of ploidy and tumor percentage values)
#' input: segmented LRR and BAF and the value for gamma_param
#' @noRd
create_distance_matrix <- function(
  s,
  dist_choice,
  gamma_param,
  uninformative_baf_threshold = 0.51,
  min_rho = 0.1,
  max_rho = 1,
  min_psi = 1,
  max_psi = 5.4
) {
  psi_pos <- seq(min_psi, max_psi, 0.05)
  rho_pos <- seq(min_rho, max_rho, 0.01)
  d <- matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) <- psi_pos
  colnames(d) <- rho_pos
  for (i in seq_along(psi_pos)) {
    psi <- psi_pos[i]
    for (j in seq_along(rho_pos)) {
      rho <- rho_pos[j]
      distance_info <- calc_distance(s,
        dist_choice, rho, psi, gamma_param,
        uninformative_baf_threshold = uninformative_baf_threshold
      )
    }
  }

  minimise <- distance_info$minimise
  return(list(distance_matrix = d, minimise = minimise))
}

#' Helper function to create the clonal distance matrix for a range of
#' rho and psi values
#' @noRd
create_distance_matrix_clonal <- function(
  segs,
  dist_choice,
  gamma_param,
  read_depth,
  siglevel_BAF,
  maxdist_BAF,
  siglevel_LogR,
  maxdist_LogR,
  uninformative_baf_threshold,
  new_bounds
) {
  psi_min <- new_bounds$psi_min
  psi_max <- new_bounds$psi_max
  rho_min <- new_bounds$rho_min
  rho_max <- new_bounds$rho_max

  s <- segs

  psi_range <- psi_max - psi_min
  rho_range <- rho_max - rho_min

  delta_psi <- psi_range / 100
  delta_rho <- rho_range / 100

  psi_pos <- seq(psi_min, psi_max, delta_psi)
  rho_pos <- seq(rho_min, rho_max, delta_rho)

  ref_seg_matrix <- matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  ref_major <- matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  ref_minor <- matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(ref_seg_matrix) <- psi_pos
  colnames(ref_seg_matrix) <- rho_pos
  rownames(ref_major) <- psi_pos
  colnames(ref_major) <- rho_pos
  rownames(ref_minor) <- psi_pos
  colnames(ref_minor) <- rho_pos

  d <- matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) <- psi_pos
  colnames(d) <- rho_pos
  for (i in seq_along(psi_pos)) {
    psi <- psi_pos[i]
    for (j in seq_along(rho_pos)) {
      rho <- rho_pos[j]

      distance_info <- calc_distance_clonal(
        s, dist_choice,
        rho, psi,
        gamma_param, read_depth,
        siglevel_BAF, maxdist_BAF,
        siglevel_LogR, maxdist_LogR,
        uninformative_baf_threshold
      )

      distance_value <- distance_info$distance_value
      max_clonal_segment <- distance_info$max_clonal_segment

      d[i, j] <- distance_value
      ref_seg_matrix[i, j] <- max_clonal_segment

      ref_major[i, j] <- distance_info$ref_maj
      ref_minor[i, j] <- distance_info$ref_min
    }
  }

  minimise <- distance_info$minimise
  return(list(
    distance_matrix = d,
    minimise = minimise,
    ref_seg_matrix = ref_seg_matrix,
    ref_major = ref_major,
    ref_minor = ref_minor
  ))
}
