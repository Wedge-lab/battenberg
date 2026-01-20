####################################################################################################
#' This function calculates a P-value, for a test where the null hypothesis is that
#' the sample was drawn from a Gaussian population with the specified mean "mu_pop".
#' @noRd
calc_Pvalue_t_twotailed <- function(
  sample_size,
  sample_mean,
  sample_SD,
  mu_pop,
  max_dist
) {
  tvar <- (sample_mean - mu_pop) * sqrt(sample_size) / sample_SD

  # We use abs(tvar) to always get the upper tail, then multiply by 2
  pval <- 2 * stats::pt(abs(tvar), df = sample_size - 1, lower.tail = FALSE)
  pval[abs(sample_mean - mu_pop) < max_dist] <- 1
  return(pval)
}

####################################################################################################
#' Helper function that calculates a binomial probability
#' @noRd
calc_binomial_prob <- function(sample_proportion, sample_size, pop_proportion) {
  p <- pmax(0, pmin(1, pop_proportion))
  x <- round(sample_proportion * sample_size)
  x <- pmax(0, pmin(sample_size, x))

  return(stats::dbinom(x, size = sample_size, prob = p))
}

####################################################################################################
#' This function calculates a log likelihood ratio where the two hypotheses are that
#' the tumour genome segment in question is "clonal".
#' The first hypothesis is the "best fit" model we can find.
#' The second hypothesis is the "second best fit" model we can find.
#' @noRd
calc_ln_likelihood_ratio <- function(LogR, BAF_req, BAF_length, BAF_size, BAF_mean, read_depth, rho, psi, gamma_param, maxdist_BAF) {
  pooled_BAF_size <- read_depth * BAF_size

  # if we don't have a value for LogR, fill in 0
  if (is.na(LogR)) {
    LogR <- 0
  }
  nMajor <- (rho - 1 + BAF_req * psi * 2^(LogR / gamma_param)) / rho
  nMinor <- (rho - 1 + (1 - BAF_req) * psi * 2^(LogR / gamma_param)) / rho


  # DCW - increase nMajor and nMinor together, to avoid impossible combinations (with negative subclonal fractions)
  if (nMinor < 0 || is.na(nMinor)) {
    if (BAF_req == 1) {
      # avoid calling infinite copy number
      nMajor <- 1000
    } else {
      nMajor <- nMajor + BAF_req * (0.01 - nMinor) / (1 - BAF_req)
      if (nMajor < 0) nMajor <- 1000
    }
    nMinor <- 0.01
  }

  if (!is.finite(nMajor)) {
    nMajor <- 0.01
  }

  # Check if there is a viable solution
  if (!is.na(BAF_req)) {
    nearest_edge <- prioritizeCopyNumbers(
      rho = rho,
      psi = psi,
      BAF_req = BAF_req,
      nMajor = nMajor,
      nMinor = nMinor,
      full = FALSE
    )
    nMaj <- nearest_edge$nMaj
    nMin <- nearest_edge$nMin
    BAF_levels <- (1 - rho + rho * nMaj) / (2 - 2 * rho + rho * (nMaj + nMin))
    index_vect <- which(is.finite(BAF_levels))
    BAF_levels <- BAF_levels[index_vect]

    if (length(BAF_levels) > 1) {
      likelihood_vect <- sapply(BAF_levels, function(x) {
        calc_binomial_prob(BAF_mean, pooled_BAF_size, x)
      })
      likelihood_vect <- sort(likelihood_vect, decreasing = TRUE)

      if ((likelihood_vect[1] > 0) && (likelihood_vect[2] > 0)) {
        ln_lratio <- log(likelihood_vect[1]) - log(likelihood_vect[2])
      } else {
        ln_lratio <- 0
      }
    } else {
      ln_lratio <- 0
    }
  } else {
    ln_lratio <- 0
  }

  return(ln_lratio)
}


#' Helper function to estimate rho from a given copy number state and it's BAF. The LogR is not used.
#' @noRd
estimate_rho <- function(LogR_value, BAF_req_value, nA_value, nB_value) {
  rho_value <- (2 * BAF_req_value - 1) / (2 * BAF_req_value - BAF_req_value * (nA_value + nB_value) - 1 + nA_value)
  return(rho_value)
}

####################################################################################################
#' Helper function to calculate psi from a copy number fit, BAF, LogR, rho and a platform gamma
#' @noRd
estimate_psi <- function(LogR_value, BAF_req_value, nA_value, nB_value, rho_value, gamma_param) {
  temp_value <- 2^(-LogR_value / gamma_param)
  temp_value <- temp_value * (2 + (rho_value * (nA_value + nB_value - 2)))
  # DCW this returns psi rather than psi_t, i.e. the average ploidy of normal and tumour cells
  temp_value <- temp_value - (2 * (1 - rho_value))
  psi_value <- temp_value / rho_value
  return(psi_value)
}


#' Function that calculates rho and psi from a given reference segment, defined by ref_seg, with copy number state nA_ref and nB_ref
#' @noRd
get_psi_rho_from_ref_seg <- function(ref_seg, s, nA_ref, nB_ref, gamma_param = 1) {
  BAF_req <- s[ref_seg, "b"]
  LogR <- s[ref_seg, "r"]

  rho <- estimate_rho(LogR, BAF_req, nA_ref, nB_ref)
  psi <- estimate_psi(LogR, BAF_req, nA_ref, nB_ref, rho, gamma_param)

  # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
  nA <- (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"] / gamma_param) * ((1 - rho) * 2 + rho * psi)) / rho
  nB <- (rho - 1 + s[, "b"] * 2^(s[, "r"] / gamma_param) * ((1 - rho) * 2 + rho * psi)) / rho
  ploidy <- sum((nA + nB) * s[, "length"]) / sum(s[, "length"])

  # TODO DEBUG
  if (rho > 0) {
    ref_segment_info <- list(psi = psi, rho = rho, ploidy = ploidy)
  } else {
    ref_segment_info <- list(psi = NA, rho = NA, ploidy = NA)
  }

  return(ref_segment_info)
}


####################################################################################################
#' This function calculates a t variate.
#' @noRd
calc_standardised_error <- function(
  LogR,
  BAF_req,
  BAF_length,
  BAF_size,
  BAF_mean,
  BAF_sd,
  rho,
  psi,
  gamma_param,
  maxdist_BAF
) {
  # if we don't have a value for LogR, fill in 0
  LogR <- ifelse(is.na(LogR), 0, LogR)

  # Pre-calculating the shared power term for clarity
  scale_factor <- psi * 2^(LogR / gamma_param)
  nMajor <- (rho - 1 + BAF_req * scale_factor) / rho
  nMinor <- (rho - 1 + (1 - BAF_req) * scale_factor) / rho

  # Floor at 0.01 (enforce "positive square")
  nMajor <- pmax(0.01, nMajor)
  nMinor <- pmax(0.01, nMinor)

  # We test 4 rounding combinations to see which matches BAF_req best
  nMaj_opts <- list(floor(nMajor), ceiling(nMajor), floor(nMajor), ceiling(nMajor))
  nMin_opts <- list(ceiling(nMinor), ceiling(nMinor), floor(nMinor), floor(nMinor))

  # Compute BAF levels for all 4 options (Vectorized)
  BAF_levels <- lapply(1:4, function(k) {
    denom <- (2 - 2 * rho + rho * (nMaj_opts[[k]] + nMin_opts[[k]]))
    (1 - rho + rho * nMaj_opts[[k]]) / pmax(denom, 1e-10)
  })

  # Find the best option for each segment (Vectorized)
  # diffs will be N x 4 matrix
  diffs <- cbind(
    abs(BAF_levels[[1]] - BAF_req),
    abs(BAF_levels[[2]] - BAF_req),
    abs(BAF_levels[[3]] - BAF_req),
    abs(BAF_levels[[4]] - BAF_req)
  )

  # Tie-breaking logic (Handle the 0.5 case for each segment)
  # This is usually for balanced regions.
  # We use max.col to find the index of the minimum difference
  best_idx <- max.col(-diffs, ties.method = "first")

  # Extract the best mu values
  # mu <- rep(0, length(BAF_req))
  # for(k in 1:4) mu[best_idx == k] <- BAF_levels[[k]][best_idx == k]
  # More R-idiomatic way:
  mu <- vapply(seq_along(best_idx), function(i) BAF_levels[[best_idx[i]]][i], numeric(1))

  # Final t-variable calculation
  is_valid <- (BAF_size > 0 & BAF_sd != 0)
  tvar <- ifelse(is_valid, (BAF_mean - mu) * sqrt(BAF_size) / BAF_sd, 0)

  return(list(included_segment = as.numeric(is_valid), tvar = tvar))
}


#' Recalculate psi_t based on rho and the available data
#'
#' @param psi A psi estimate
#' @param rho A rho estimate
#' @param platform_gamma The platform specific LogR scaling parameter
#' @param lrrsegmented Segmented LogR, a vector with just the values
#' @param segBAF_table Segmented BAF, the full table
#' @param siglevel_BAF Significance level when testing wether a segment is clonal or subclonal given a rho/psi combination, parameter is used in \code{is_segment_clonal}
#' @param maxdist_BAF Max distance BAF is allowed to be away from the copy number solution before we don't trust the value and overrule a p-value, parameter required when determining the clonal status of a segment in \code{is_segment_clonal}
#' @param include_subcl_segments Boolean flag, supply TRUE if subclonal segments should be included when calculating psi_t, supply FALSE if only clonal segments should be included (default: TRUE)
#' @noRd
recalc_psi_t <- function(psi, rho, gamma_param, lrrsegmented, segBAF_table, siglevel_BAF, maxdist_BAF, include_subcl_segments = TRUE) {
  # Create segments of constant BAF/LogR
  s <- get_segment_info(lrrsegmented[rownames(segBAF_table)], segBAF_table)
  # Make sure no segment of length 1 remains - TODO: this should not occur and needs to be prevented upstream
  s <- s[s[, 3] > 1, ]

  # Check which segments are clonal with this rho/psi configuration
  segment_info <- is_segment_clonal(
    LogR = s[, "r"],
    BAF_req = s[, "b"],
    BAF_length = s[, "length"],
    BAF_size = s[, "size"],
    BAF_mean = s[, "mean"],
    BAF_sd = s[, "sd"],
    rho = rho,
    psi = psi,
    gamma_param = gamma_param,
    siglevel_BAF = siglevel_BAF,
    maxdist_BAF = maxdist_BAF
  )

  # Include this segment if we want to include all segments,
  # or if we don't want subclonal segments include it only if its clonal
  keep_mask <- if (include_subcl_segments) rep(TRUE, nrow(s)) else segment_info$is_clonal

  if (!any(keep_mask)) {
    return(NA)
  }

  nMaj <- segment_info$nMaj[keep_mask]
  nMin <- segment_info$nMin[keep_mask]
  s_r <- s[keep_mask, "r"]
  s_len <- s[keep_mask, "length"]

  # Calculate psi_t for each segment and then the weighted average
  psi_t_vec <- calc_psi_t(nMaj + nMin, s_r, rho, gamma_param)
  psi_t <- collapse::fsum(psi_t_vec * s_len) / collapse::fsum(s_len)

  return(psi_t)
}

#' Calculate psi based on a reference segment and its associated logr
#'
#' @param total_cn Integer representing the total clonal copynumber (i.e. nMajor+nMinor)
#' @param r The LogR of the segment with the total_cn copy number
#' @param rho A cellularity estimate
#' @param gamma_param Platform gamma parameter
#' @author sd11
#' @export
calc_psi_t <- function(total_cn, r, rho, gamma_param) {
  psi <- (rho * (total_cn) + 2 - 2 * rho) / (2^(r / gamma_param))
  psi_t <- (psi - 2 * (1 - rho)) / rho
  return(psi_t)
}


# Optimized Batch version of the t-test logic
calc_batch_standardised_errors <- function(s, rho, psi, gamma_param) {
  # s contains columns: r (LogR), b (BAF_req), length, size, mean, sd

  scale <- psi * 2^(s[, "r"] / gamma_param)
  nMajor_raw <- (rho - 1 + s[, "b"] * scale) / rho
  nMinor_raw <- (rho - 1 + (1 - s[, "b"]) * scale) / rho

  # Vectorized floor at 0.01
  nMajor <- pmax(0.01, nMajor_raw)
  nMinor <- pmax(0.01, nMinor_raw)

  # Instead of a 4-item list per segment, we do 4 separate vector calculations
  # This is where the massive speedup happens
  nMaj_opts <- list(floor(nMajor), ceiling(nMajor), floor(nMajor), ceiling(nMajor))
  nMin_opts <- list(ceiling(nMinor), ceiling(nMinor), floor(nMinor), floor(nMinor))

  # Calculate BAF levels for all 4 possibilities across all segments simultaneously
  BAF_levels <- lapply(1:4, function(k) {
    denom <- (2 - 2 * rho + rho * (nMaj_opts[[k]] + nMin_opts[[k]]))
    (1 - rho + rho * nMaj_opts[[k]]) / denom
  })

  # Vectorized "which.min(abs(BAF_levels - BAF_req))"
  # We find the distance for all 4 options
  diffs <- cbind(
    abs(BAF_levels[[1]] - s[, "b"]),
    abs(BAF_levels[[2]] - s[, "b"]),
    abs(BAF_levels[[3]] - s[, "b"]),
    abs(BAF_levels[[4]] - s[, "b"])
  )

  # Pick the best index for every segment at once
  best_idx <- max.col(-diffs) # max of negative is min

  # Map the best mu values
  mu <- mapply(function(row, col) BAF_levels[[col]][row], seq_len(nrow(s)), best_idx)

  # Final t-variable calculation
  is_valid <- s[, "size"] > 0 & s[, "sd"] != 0
  tvar <- ifelse(is_valid, (s[, "mean"] - mu) * sqrt(s[, "size"]) / s[, "sd"], 0)

  return(tvar)
}
