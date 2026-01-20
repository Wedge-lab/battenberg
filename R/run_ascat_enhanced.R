#' Key optimizations:
#' 1. Early termination after first good solution (like original)
#' 2. Vectorized distance calculations
#' 3. Optimized constraint checking
#' 4. Smart search ordering (best regions first)
#' 5. Reduced memory allocations
runASCAT_enhanced <- function(
  lrr, baf, lrrsegmented, bafsegmented, chromosomes, dist_choice,
  distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA,
  cnaStatusFile = "copynumber_solution_status.txt", gamma = 0.55,
  allow100percent, reliabilityFile = NA, min_ploidy = 1.6, max_ploidy = 4.8,
  min_rho = 0.1, max_rho = 1.0, min_goodness = 63,
  uninformative_baf_threshold = 0.51, chr_names, analysis = "paired",
  smart_ordering = TRUE, early_termination = TRUE, verbose = TRUE, nthreads = 1
) {
  start_time <- Sys.time()

  # 1. Setup Data Processing
  ch <- chromosomes
  b <- bafsegmented
  r <- lrrsegmented[names(bafsegmented)]

  dist_min_psi <- max(min_ploidy - 0.6, 0)
  dist_max_psi <- max_ploidy + 0.6
  dist_min_rho <- max(min_rho - 0.03, 0.05)
  dist_max_rho <- max_rho + 0.03

  # 2. Create Segments & Distance Matrix
  s <- make_segments(r, b)
  dist_matrix_info <- create_distance_matrix(s, dist_choice, gamma,
    uninformative_baf_threshold = uninformative_baf_threshold,
    min_psi = dist_min_psi, max_psi = dist_max_psi,
    min_rho = dist_min_rho, max_rho = dist_max_rho,
    nthreads = nthreads
  )
  d <- dist_matrix_info$distance_matrix

  # Theoretical maximum distance (weighted by length)
  TheoretMaxdist <- collapse::fsum(rep(0.25, nrow(s)) * s[, "length"],
    na.rm = TRUE
  )

  minimise <- dist_matrix_info$minimise

  log_debug("--- Debug: Grid and Segments ---")
  log_debug("Number of segments created: {nrow(s)}")
  log_debug("Distance matrix dimensions: {nrow(d)} x: {ncol(d)}")
  log_debug("Theoretical Max Distance: {round(TheoretMaxdist, 4)}")
  if (!minimise) d <- -d

  # 3. Pre-compute Search Parameters
  rho_values <- as.numeric(colnames(d))
  psi_values <- as.numeric(rownames(d))
  s_length <- s[, "length"]
  s_b <- s[, "b"]
  s_r <- s[, "r"]
  total_length <- collapse::fsum(s_length)

  # Pre-compute masks for calculate_solution_fast
  baf_mask <- s_b != 0.5
  denom_abb <- collapse::fsum(s_length[baf_mask])

  # Get search matrix (i, j)
  search_order <- create_smart_search_order(d, smart_ordering, verbose)
  total_points_in_grid <- nrow(search_order)

  # 4. Main Search Loop
  nropt <- 0
  optima <- list()
  localmin_vals <- numeric()
  points_checked <- 0

  if (total_points_in_grid > 0) {
    for (idx in seq_len(total_points_in_grid)) {
      i <- search_order[idx, 1]
      j <- search_order[idx, 2]
      m <- d[i, j]
      points_checked <- points_checked + 1

      if (is_local_minimum_fast(d, i, j, m)) {
        solution <- calculate_solution_fast(
          psi_values[i], rho_values[j], s_b, s_r, s_length, total_length,
          gamma, min_ploidy, max_ploidy, min_rho, max_rho,
          min_goodness, m, TheoretMaxdist, minimise, allow100percent,
          baf_mask = baf_mask, denom_abb = denom_abb
        )

        if (!solution_is_null(solution)) {
          nropt <- nropt + 1
          # Store as vector for consistency with original optima extraction
          optima[[nropt]] <- c(m, i, j, solution$ploidy, solution$goodness)
          localmin_vals[nropt] <- m

          if (verbose) {
            log_info("Found solution {nropt} at point {points_checked}: rho={round(rho_values[j], 3)}, psi={round(psi_values[i], 3)}")
          }

          if (early_termination && solution$goodness >= (min_goodness + 5)) break
        }
      }
      if (verbose && points_checked %% 5000 == 0) log_info("Progress: {points_checked} points checked")
    }
  }

  # 5. Handle 100% Aberrant Fallback
  if (allow100percent && nropt == 0) {
    if (verbose) log_info("Trying 100% aberrant solutions...")
    d_mod <- d
    d_mod[, rho_values <= 1] <- 1e20
    search_order_100 <- create_smart_search_order(d_mod, smart_ordering, FALSE)

    if (nrow(search_order_100) > 0) {
      for (idx in seq_len(nrow(search_order_100))) {
        i <- search_order_100[idx, 1]
        j <- search_order_100[idx, 2]
        m <- d_mod[i, j]
        if (is_local_minimum_fast(d_mod, i, j, m)) {
          solution <- calculate_solution_fast(
            psi_values[i], rho_values[j], s_b, s_r, s_length, total_length, gamma,
            min_ploidy, max_ploidy, min_rho, max_rho,
            min_goodness, m, TheoretMaxdist, minimise, allow100percent,
            baf_mask = baf_mask, denom_abb = denom_abb,
            skip_zero_check = FALSE
          )
          if (!solution_is_null(solution)) {
            nropt <- 1
            optima[[1]] <- c(m, i, j, solution$ploidy, solution$goodness)
            localmin_vals[1] <- m
            break
          }
        }
      }
    }
  }

  optimization_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # 6. Select Best Solution & Collect Sunrise Plot Data
  if (nropt > 0) {
    data.table::fwrite(list(paste0(nropt, " copy number solutions found")), cnaStatusFile)

    optlim <- sort(localmin_vals)[1]
    psi_opt1_plot <- numeric()
    rho_opt1_plot <- numeric()

    # Original logic: collect all solutions that share the global minimum distance
    for (idx in seq_along(optima)) {
      if (optima[[idx]][1] == optlim) {
        psi_opt1 <- psi_values[optima[[idx]][2]]
        rho_opt1 <- min(rho_values[optima[[idx]][3]], 1.0)
        ploidy_opt1 <- optima[[idx]][4]
        goodness_of_fit_opt1 <- optima[[idx]][5]

        psi_opt1_plot <- c(psi_opt1_plot, psi_opt1)
        rho_opt1_plot <- c(rho_opt1_plot, rho_opt1)
      }
    }
  } else {
    data.table::fwrite(list("no copy number solutions found"), cnaStatusFile)
    return(list(
      psi = NA, rho = NA, ploidy = NA,
      convergence_info = list(
        converged = FALSE, n_solutions_found = 0,
        optimization_time = optimization_time, points_checked = points_checked,
        search_efficiency = points_checked / total_points_in_grid
      )
    ))
  }

  # Use the extracted "best" values for the final vectors
  rho <- rho_opt1
  psi <- psi_opt1
  ploidy <- ploidy_opt1
  goodness_of_fit <- goodness_of_fit_opt1

  # 7. Final Back-transformation
  mult <- 2^(r / gamma) * ((1 - rho) * 2 + rho * psi)
  nAfull <- (rho - 1 - (b - 1) * mult) / rho
  nBfull <- (rho - 1 + b * mult) / rho
  nA <- pmax(round(nAfull), 0)
  nB <- pmax(round(nBfull), 0)

  rBacktransform <- gamma * log(
    (rho * (nA + nB) +
      (1 - rho) * 2) / ((1 - rho) * 2 + rho * psi),
    2
  )
  bBacktransform <- (1 - rho + rho * nB) / (2 - 2 * rho + rho * (nA + nB))

  # Logic check: ensures reliability metrics are identical to original source
  # Logic check: ensures reliability metrics are identical to original source
  rDiff <- 1 - abs(rBacktransform - r) / abs(r)
  rConf <- ifelse(abs(rBacktransform) > 0.15,
    pmin(100, pmax(0, 100 * rDiff)), NA
  )
  bDiff <- 1 - abs(bBacktransform - b) / abs(b - 0.5)
  bConf <- ifelse(bBacktransform != 0.5,
    pmin(100, pmax(0, ifelse(b == 0.5, 100, 100 * bDiff))), NA
  )

  if (!is.na(reliabilityFile)) {
    data.table::fwrite(
      data.frame(
        segmentedBAF = b, backTransformedBAF = bBacktransform,
        confidenceBAF = bConf, segmentedR = r,
        backTransformedR = rBacktransform, confidenceR = rConf,
        nA = nA, nB = nB, nAfull = nAfull, nBfull = nBfull
      ),
      reliabilityFile,
      sep = ",", row.names = FALSE
    )
  }

  # 8. Plotting
  if (analysis == "paired" && !is.na(distancepng)) {
    grDevices::png(filename = distancepng, width = 1000, height = 1000, res = 150, type = "cairo")
    ASCAT::ascat.plotSunrise(-d, psi_opt1_plot, rho_opt1_plot, minimise)
    grDevices::dev.off()
  }

  if (!is.na(copynumberprofilespng)) {
    grDevices::png(filename = copynumberprofilespng, width = 2000, height = 500, res = 200, type = "cairo")
    ASCAT::ascat.plotAscatProfile(
      n1all = nA, n2all = nB, heteroprobes = TRUE, ploidy = ploidy,
      rho = rho, goodness_of_fit = goodness_of_fit, nonaberrant = FALSE,
      ch = ch, lrr = lrr, bafsegmented = bafsegmented, chrs = chr_names
    )
    grDevices::dev.off()
  }

  if (!is.na(nonroundedprofilepng)) {
    grDevices::png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200, type = "cairo")
    ASCAT::ascat.plotNonRounded(
      ploidy = ploidy, rho = rho, goodness_of_fit = goodness_of_fit,
      nonaberrant = FALSE, nAfull = nAfull, nBfull = nBfull,
      bafsegmented = bafsegmented, ch = ch, lrr = lrr, chrs = chr_names
    )
    grDevices::dev.off()
  }

  return(list(
    psi = psi, rho = rho, ploidy = ploidy,
    convergence_info = list(
      converged = TRUE,
      n_solutions_found = nropt,
      optimization_time = optimization_time,
      points_checked = points_checked,
      search_efficiency = points_checked / total_points_in_grid
    )
  ))
}

create_smart_search_order <- function(d, smart_ordering, verbose) {
  idx_mat <- which(is.finite(d), arr.ind = TRUE)
  if (nrow(idx_mat) == 0) {
    return(matrix(0, 0, 2))
  }

  nr <- nrow(d)
  nc <- ncol(d)
  # Original Battenberg border logic: 4:(nr-3)
  # We only apply it if the matrix is large enough to have an interior
  if (nr >= 7 && nc >= 7) {
    keep <- idx_mat[, 1] >= 4 & idx_mat[, 1] <= (nr - 3) &
      idx_mat[, 2] >= 4 & idx_mat[, 2] <= (nc - 3)
    # If the border filter leaves points, use them; otherwise keep original (edge case)
    if (any(keep)) idx_mat <- idx_mat[keep, , drop = FALSE]
  }

  if (smart_ordering) {
    # Extract distances via matrix indexing (no loop)
    distances <- d[idx_mat]
    idx_mat <- idx_mat[order(distances), ]
  }
  return(idx_mat)
}

#' Fast solution calculation (vectorized and optimized)
calculate_solution_fast <- function(
  psi, rho, s_b, s_r, s_length, total_length, gamma,
  min_ploidy, max_ploidy, min_rho, max_rho,
  min_goodness, distance_value, TheoretMaxdist, minimise,
  allow100percent, baf_mask, denom_abb, skip_zero_check = FALSE
) {
  # Guard against rho = 0 to prevent Inf
  safe_rho <- pmax(rho, 1e-6)

  # Constraint pre-check
  if (psi < min_ploidy || psi > max_ploidy || rho < min_rho || rho > max_rho) {
    return(NULL)
  }

  # Vectorized calculation
  multiplier <- 2^(s_r / gamma) * ((1 - safe_rho) * 2 + safe_rho * psi)
  nA <- (safe_rho - 1 - (s_b - 1) * multiplier) / safe_rho
  nB <- (safe_rho - 1 + s_b * multiplier) / safe_rho

  # Ploidy check
  ploidy <- collapse::fsum((nA + nB) * s_length) / total_length
  if (is.na(ploidy) || ploidy < min_ploidy || ploidy > max_ploidy) {
    return(NULL)
  }

  # Goodness check
  goodness_of_fit <- if (minimise) {
    (1 - distance_value / TheoretMaxdist) * 100
  } else {
    -distance_value / TheoretMaxdist * 100
  }
  if (is.na(goodness_of_fit) || goodness_of_fit < min_goodness) {
    return(NULL)
  }

  if (!skip_zero_check && !allow100percent) {
    nA_r <- round(nA)
    nB_r <- round(nB)
    # Edge case: sum(s_length[logical]) can be 0 if no indices match
    percentzero <- (collapse::fsum(s_length[which(nA_r == 0)]) +
      collapse::fsum(s_length[which(nB_r == 0)])) / total_length

    perczeroAbb <- 0
    if (denom_abb > 0) {
      # Use which() to avoid NA issues in logical indexing
      # Use which() to avoid NA issues in logical indexing
      perczeroAbb <- (collapse::fsum(s_length[which(baf_mask & nA_r == 0)]) +
        collapse::fsum(s_length[which(baf_mask & nB_r == 0)])) /
        denom_abb
    }
    if (!(percentzero > 0.01 || perczeroAbb > 0.1)) {
      return(NULL)
    }
  }

  return(list(psi = psi, rho = min(rho, 1.0), ploidy = ploidy, goodness = goodness_of_fit))
}

#' Fast local minimum check (optimized version of original 7x7)
is_local_minimum_fast <- function(d, i, j, center_value) {
  # Check 7x7 neighborhood (same as original)
  i_min <- i - 3
  i_max <- i + 3
  j_min <- j - 3
  j_max <- j + 3

  # Bounds checking
  if (i_min < 1 || i_max > nrow(d) || j_min < 1 || j_max > ncol(d)) {
    return(FALSE)
  }

  # Extract neighborhood
  neighborhood <- d[i_min:i_max, j_min:j_max]

  # Set center to maximum to exclude it from minimum check
  neighborhood[4, 4] <- max(neighborhood, na.rm = TRUE)

  # Check if center is local minimum
  return(min(neighborhood, na.rm = TRUE) > center_value)
}

solution_is_null <- function(sol) {
  return(is.null(sol) || is.na(sol$ploidy) || is.na(sol$goodness))
}
