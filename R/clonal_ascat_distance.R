#' Prepare segment data for optimization
#'
#' Filters informative segments and extracts vectors
#' @noRd
prepare_clonal_segments <- function(segs, uninformative_baf_threshold) {
  # segs has columns: r, b, length, size, mean, sd
  informative_idx <- segs[, "b"] > uninformative_baf_threshold
  if (!any(informative_idx)) {
    return(NULL)
  }

  list(
    r = segs[informative_idx, "r"],
    b = segs[informative_idx, "b"],
    len = segs[informative_idx, "length"],
    size = segs[informative_idx, "size"],
    mean = segs[informative_idx, "mean"],
    sd = segs[informative_idx, "sd"],
    genome_size = sum(segs[informative_idx, "size"]),
    total_segs = length(segs[informative_idx, "b"])
  )
}

#' Numeric-only version of calc_distance_clonal
#' Uses pre-extracted vectors for speed
#' @noRd
calc_clonal_distance_numeric <- function(
  seg_list, dist_choice, rho, psi, gamma_param, read_depth,
  siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR
) {
  if (is.null(seg_list)) {
    return(list(
      distance_value = 0,
      minimise = FALSE,
      max_clonal_segment = 0, # Return 0 index if no informative segments
      ref_maj = NA,
      ref_min = NA
    ))
  }

  # Call vectorized is_segment_clonal with vectors from the list
  seg_info <- is_segment_clonal(
    LogR = seg_list$r, BAF_req = seg_list$b, BAF_length = seg_list$len,
    BAF_size = seg_list$size, BAF_mean = seg_list$mean, BAF_sd = seg_list$sd,
    rho = rho, psi = psi, gamma_param = gamma_param,
    siglevel_BAF = siglevel_BAF, maxdist_BAF = maxdist_BAF
  )

  # Calculate Standardised Error (Vectorized)
  err_info <- calc_standardised_error(
    LogR = seg_list$r, BAF_req = seg_list$b, BAF_length = seg_list$len,
    BAF_size = seg_list$size, BAF_mean = seg_list$mean, BAF_sd = seg_list$sd,
    rho = rho, psi = psi, gamma_param = gamma_param,
    maxdist_BAF = maxdist_BAF
  )

  # Extract results
  is_clonal <- seg_info$is_clonal
  is_balanced <- seg_info$balanced
  nMaj <- seg_info$nMaj
  nMin <- seg_info$nMin
  tvar_sq <- err_info$tvar^2

  # Max Clonal Segment logic
  # Determine best match index (relative to the SUBSETTED list)
  potential_indices <- which(is_clonal & !is_balanced)

  max_idx_local <- 0
  ref_maj <- NA
  ref_min <- NA

  if (length(potential_indices) > 0) {
    # Find index in the subset
    local_best <- potential_indices[which.max(seg_list$size[potential_indices])]
    max_idx_local <- local_best # This is the index in the *informative* subset
    ref_maj <- nMaj[local_best]
    ref_min <- nMin[local_best]
  }

  # Compute Distance
  dist_val <- 0
  minimise <- FALSE
  dc <- as.character(dist_choice)

  if (dc == "0") {
    clonal_genome_size <- sum(seg_list$size[is_clonal])
    dist_val <- clonal_genome_size / pmax(seg_list$genome_size, 1e-10)
    minimise <- FALSE
  } else if (dc == "1") {
    n_inc <- sum(err_info$included_segment)
    dist_val <- sum(tvar_sq) / pmax(n_inc, 1)
    minimise <- TRUE
  } else if (dc == "2") {
    b_diff_sq <- (seg_list$b - seg_list$mean)^2
    dist_val <- sum(b_diff_sq) / pmax(seg_list$total_segs, 1)
    minimise <- TRUE
  } else if (dc == "3") {
    b_diff_sq <- (seg_list$b - seg_list$mean)^2
    dist_val <- sum(seg_list$size * b_diff_sq) / pmax(seg_list$genome_size, 1e-10)
    minimise <- TRUE
  } else if (dc == "4") {
    ln_lratio <- calc_ln_likelihood_ratio(
      LogR = seg_list$r, BAF_req = seg_list$b, BAF_length = seg_list$len,
      BAF_size = seg_list$size, BAF_mean = seg_list$mean, read_depth = read_depth,
      rho = rho, psi = psi, gamma_param = gamma_param,
      maxdist_BAF = maxdist_BAF
    )
    dist_val <- sum(ln_lratio)
    minimise <- FALSE
  }

  return(list(
    distance_value     = dist_val,
    minimise           = minimise,
    max_clonal_segment = max_idx_local, # Note: this is local index!
    ref_maj            = ref_maj,
    ref_min            = ref_min
  ))
}

# Kept for backward compatibility if needed, but unused in optimized path
calc_distance_clonal <- function(
  segs, dist_choice, rho, psi, gamma_param, read_depth,
  siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR,
  uninformative_baf_threshold
) {
  # Wrap the new logic: prepare then calc
  seg_list <- prepare_clonal_segments(segs, uninformative_baf_threshold)

  # Need to map local index back to global index for this legacy wrapper function
  res <- calc_clonal_distance_numeric(
    seg_list, dist_choice, rho, psi, gamma_param, read_depth,
    siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR
  )

  # Map back the index if valid
  if (res$max_clonal_segment > 0) {
    informative_idx <- which(segs[, "b"] > uninformative_baf_threshold)
    res$max_clonal_segment <- informative_idx[res$max_clonal_segment]
  }

  return(res)
}

####################################################################################################
#' function to create the distance matrix (distance for a range of ploidy and tumor percentage values)
#' input: segmented LRR and BAF and the value for gamma_param
#' @noRd
create_distance_matrix <- function(
  s, dist_choice, gamma_param, uninformative_baf_threshold = 0.51,
  min_rho = 0.1, max_rho = 1, min_psi = 1, max_psi = 5.4, nthreads = 1
) {
  psi_pos <- seq(min_psi, max_psi, 0.05)
  rho_pos <- seq(min_rho, max_rho, 0.01)

  log_info("DEBUG: create_distance_matrix called with nthreads={nthreads}")
  log_info("DEBUG: Grid size: {length(psi_pos)}x{length(rho_pos)} \\
           ({length(psi_pos)*length(rho_pos)} iterations)")


  # PRE-EXTRACT COLUMNS (Massive speedup: stop looking up "s[,col]" inside loops)
  s_r <- s[, "r"]
  s_b <- s[, "b"]
  s_len <- s[, "length"]
  s_size <- s[, "size"]
  s_mean <- s[, "mean"]
  s_sd <- s[, "sd"]

  logR_term <- 2^(s_r / gamma_param)

  # Define the row calculation function
  calc_row <- function(psi) {
    scale_factor <- psi * logR_term
    vapply(rho_pos, function(rho) {
      nMaj_raw <- (rho - 1 + s_b * scale_factor) / rho
      nMin_raw <- (rho - 1 + (1 - s_b) * scale_factor) / rho
      nM_J <- pmax(0.01, nMaj_raw)
      nM_N <- pmax(0.01, nMin_raw)

      nMaj_opts <- list(floor(nM_J), ceiling(nM_J), floor(nM_J), ceiling(nM_J))
      nMin_opts <- list(ceiling(nM_N), ceiling(nM_N), floor(nM_N), floor(nM_N))

      best_dist <- rep(Inf, length(s_b))
      best_mu <- rep(0, length(s_b))

      for (k in 1:4) {
        denom <- (2 - 2 * rho + rho * (nMaj_opts[[k]] + nMin_opts[[k]]))
        mu_opt <- (1 - rho + rho * nMaj_opts[[k]]) / pmax(denom, 1e-10)
        dist_to_b <- abs(mu_opt - s_b)
        better <- !is.na(dist_to_b) & dist_to_b < best_dist
        best_dist[better] <- dist_to_b[better]
        best_mu[better] <- mu_opt[better]
      }

      is_valid <- s_size > 0 & s_sd != 0
      tvar <- ifelse(is_valid, (s_mean - best_mu) * sqrt(s_size) / s_sd, 0)
      return(collapse::fsum(tvar^2 * s_len))
    }, FUN.VALUE = numeric(1))
  }

  if (nthreads > 1) {
    # Parallel execution
    rows <- parallel::mclapply(psi_pos, calc_row, mc.cores = nthreads)
    d <- do.call(rbind, rows)
  } else {
    # Serial execution
    d <- matrix(nrow = length(psi_pos), ncol = length(rho_pos))
    for (i in seq_along(psi_pos)) {
      d[i, ] <- calc_row(psi_pos[i])
    }
  }

  rownames(d) <- psi_pos
  colnames(d) <- rho_pos

  return(list(distance_matrix = d, minimise = TRUE))
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
  new_bounds,
  nthreads = 1
) {
  psi_min <- new_bounds$psi_min
  psi_max <- new_bounds$psi_max
  rho_min <- new_bounds$rho_min
  rho_max <- new_bounds$rho_max

  psi_range <- psi_max - psi_min
  rho_range <- rho_max - rho_min

  delta_psi <- psi_range / 100
  delta_rho <- rho_range / 100

  psi_pos <- seq(psi_min, psi_max, delta_psi)
  rho_pos <- seq(rho_min, rho_max, delta_rho)

  # Define calculation for a single psi (row)
  # Precompute segment invariants ONCE
  seg_list <- prepare_clonal_segments(segs, uninformative_baf_threshold)

  # For mapping back indices later
  informative_indices <- which(segs[, "b"] > uninformative_baf_threshold)

  # Define calculation for a single psi (row)
  calc_row <- function(psi) {
    len_rho <- length(rho_pos)
    d_row <- numeric(len_rho)
    r_seg_row <- numeric(len_rho)
    r_maj_row <- numeric(len_rho)
    r_min_row <- numeric(len_rho)

    for (j in seq_along(rho_pos)) {
      rho <- rho_pos[j]

      # Use the optimized numeric kernel
      distance_info <- calc_clonal_distance_numeric(
        seg_list, dist_choice, rho, psi, gamma_param, read_depth,
        siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR
      )

      d_row[j] <- distance_info$distance_value

      # Map local index back to global index for "max_clonal_segment"
      local_idx <- distance_info$max_clonal_segment
      if (local_idx > 0 && length(informative_indices) >= local_idx) {
        r_seg_row[j] <- informative_indices[local_idx]
      } else {
        r_seg_row[j] <- 0
      }

      r_maj_row[j] <- distance_info$ref_maj
      r_min_row[j] <- distance_info$ref_min
    }
    return(list(d = d_row, r_seg = r_seg_row, r_maj = r_maj_row, r_min = r_min_row))
  }

  # Execute
  if (nthreads > 1) {
    res_list <- parallel::mclapply(psi_pos, calc_row, mc.cores = nthreads)
  } else {
    res_list <- lapply(psi_pos, calc_row)
  }

  # Assemble matrices
  d <- do.call(rbind, lapply(res_list, `[[`, "d"))
  ref_seg_matrix <- do.call(rbind, lapply(res_list, `[[`, "r_seg"))
  ref_major <- do.call(rbind, lapply(res_list, `[[`, "r_maj"))
  ref_minor <- do.call(rbind, lapply(res_list, `[[`, "r_min"))

  rownames(d) <- psi_pos
  colnames(d) <- rho_pos
  rownames(ref_seg_matrix) <- psi_pos
  colnames(ref_seg_matrix) <- rho_pos
  rownames(ref_major) <- psi_pos
  colnames(ref_major) <- rho_pos
  rownames(ref_minor) <- psi_pos
  colnames(ref_minor) <- rho_pos

  # Determine minimise flag (constant for all iterations)
  # We can just check the first combination
  # Determine minimise flag (constant for all iterations)
  # We can just check the first combination using the old function or new kernel
  # Use new kernel for consistency
  temp_info <- calc_clonal_distance_numeric(
    seg_list, dist_choice, rho_pos[1], psi_pos[1], gamma_param, read_depth,
    siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR
  )
  minimise <- temp_info$minimise
  return(list(
    distance_matrix = d,
    minimise = minimise,
    ref_seg_matrix = ref_seg_matrix,
    ref_major = ref_major,
    ref_minor = ref_minor
  ))
}
