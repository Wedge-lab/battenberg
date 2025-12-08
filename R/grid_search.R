#' Key optimizations:
#' 1. Early termination after first good solution (like original)
#' 2. Vectorized distance calculations
#' 3. Optimized constraint checking
#' 4. Smart search ordering (best regions first)
#' 5. Reduced memory allocations
runASCAT_enhanced = function(lrr, baf, lrrsegmented, bafsegmented, chromosomes, dist_choice,
                             distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA,
                             cnaStatusFile = "copynumber_solution_status.txt", gamma = 0.55, allow100percent,
                             reliabilityFile=NA, min.ploidy=1.6, max.ploidy=4.8, min.rho=0.1, max.rho=1.0,
                             min.goodness=63, uninformative_BAF_threshold = 0.51, chr.names, analysis="paired",
                             smart_ordering = TRUE, early_termination = TRUE, verbose = TRUE) {

  start_time <- Sys.time()

  # Setup data processing (IDENTICAL to original)
  ch = chromosomes
  b = bafsegmented
  r = lrrsegmented[names(bafsegmented)]

  dist_min_psi = max(min.ploidy-0.6, 0)
  dist_max_psi = max.ploidy+0.6
  dist_min_rho = max(min.rho-0.03, 0.05)
  dist_max_rho = max.rho+0.03

  s = ASCAT::make_segments(r,b)
  dist_matrix_info <- create_distance_matrix(s, dist_choice, gamma,
                                           uninformative_BAF_threshold=uninformative_BAF_threshold,
                                           min_psi=dist_min_psi, max_psi=dist_max_psi,
                                           min_rho=dist_min_rho, max_rho=dist_max_rho)
  d = dist_matrix_info$distance_matrix
  minimise = dist_matrix_info$minimise

  TheoretMaxdist = sum(rep(0.25,dim(s)[1]) * s[,"length"],na.rm=T)

  if( !(minimise) ) {
    d = -d
  }

  if (verbose) {
    cat("Optimized Battenberg with smart ordering and early termination...\n")
  }

  # Pre-compute values for speed
  rho_values <- as.numeric(colnames(d))
  psi_values <- as.numeric(rownames(d))
  s_length <- s[,"length"]
  s_b <- s[,"b"]
  s_r <- s[,"r"]
  total_length <- sum(s_length)

  # Create search order - most promising regions first
  search_order <- create_smart_search_order(d, smart_ordering, verbose)

  # OPTIMIZED SEARCH with early termination
  nropt = 0
  localmin = NULL
  optima = list()
  points_checked = 0

  for (idx in 1:length(search_order)) {
    point <- search_order[[idx]]
    i <- point$i
    j <- point$j
    points_checked <- points_checked + 1

    m = d[i,j]
    if (!is.finite(m)) next

    # Fast local minimum check (7x7 like original for speed)
    if (is_local_minimum_fast(d, i, j, m)) {
      psi = psi_values[i]
      rho = rho_values[j]

      # Fast solution calculation
      solution <- calculate_solution_fast(psi, rho, s_b, s_r, s_length, total_length, gamma,
                                         min.ploidy, max.ploidy, min.rho, max.rho,
                                         min.goodness, m, TheoretMaxdist, minimise, allow100percent)

      if (!is.null(solution)) {
        nropt = nropt + 1
        optima[[nropt]] = c(m, i, j, solution$ploidy, solution$goodness)
        localmin[nropt] = m

        if (verbose) {
          cat("Found solution", nropt, "at point", points_checked, "/", length(search_order),
              ": rho=", round(solution$rho, 3), ", psi=", round(solution$psi, 3),
              ", goodness=", round(solution$goodness, 2), "\n")
        }

        # Early termination if we found a good solution
        if (early_termination && solution$goodness >= (min.goodness + 5)) {
          if (verbose) cat("Early termination - found high quality solution\n")
          break
        }
      }
    }

    # Progress update
    if (verbose && points_checked %% 2000 == 0) {
      cat("Progress:", points_checked, "/", length(search_order), "points checked\n")
    }
  }

  # Handle 100% aberrant case (only if no solutions found)
  if (allow100percent & nropt == 0) {
    if (verbose) cat("Trying 100% aberrant solutions...\n")

    cold = which(rho_values > 1)
    d_modified <- d
    d_modified[,cold] = 1E20

    # Use same optimized search for 100% case
    search_order_100 <- create_smart_search_order(d_modified, smart_ordering, FALSE)

    for (idx in 1:length(search_order_100)) {
      point <- search_order_100[[idx]]
      i <- point$i
      j <- point$j

      m = d_modified[i,j]
      if (!is.finite(m)) next

      if (is_local_minimum_fast(d_modified, i, j, m)) {
        psi = psi_values[i]
        rho = rho_values[j]

        solution <- calculate_solution_fast(psi, rho, s_b, s_r, s_length, total_length, gamma,
                                           min.ploidy, max.ploidy, min.rho, max.rho,
                                           min.goodness, m, TheoretMaxdist, minimise, allow100percent,
                                           skip_zero_check = TRUE)

        if (!is.null(solution)) {
          nropt = nropt + 1
          optima[[nropt]] = c(m, i, j, solution$ploidy, solution$goodness)
          localmin[nropt] = m
          break  # Early termination for 100% case too
        }
      }
    }
  }

  optimization_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Process results (IDENTICAL to original logic)
  psi_opt1_plot = vector(mode="numeric")
  rho_opt1_plot = vector(mode="numeric")

  if (nropt>0) {
    write.table(paste(nropt, " copy number solutions found", sep=""), file=cnaStatusFile, quote=F, col.names=F, row.names=F)
    optlim = sort(localmin)[1]

    for (i in 1:length(optima)) {
      if(optima[[i]][1] == optlim) {
        psi_opt1 = psi_values[optima[[i]][2]]
        rho_opt1 = rho_values[optima[[i]][3]]
        if(rho_opt1 > 1) {
          rho_opt1 = 1
        }
        ploidy_opt1 = optima[[i]][4]
        goodnessOfFit_opt1 = optima[[i]][5]
        psi_opt1_plot = c(psi_opt1_plot, psi_opt1)
        rho_opt1_plot = c(rho_opt1_plot, rho_opt1)
      }
    }
  } else {
    write.table(paste("no copy number solutions found", sep=""), file=cnaStatusFile, quote=F, col.names=F, row.names=F)
    if (verbose) cat("No suitable copy number solution found\n")
    psi = NA
    ploidy = NA
    rho = NA
    psi_opt1_plot = -1
    rho_opt1_plot = -1

    return(list(
      psi = psi,
      rho = rho,
      ploidy = ploidy,
      convergence_info = list(
        converged = FALSE,
        optimization_time = optimization_time,
        points_checked = points_checked
      )
    ))
  }

  if (verbose) {
    cat("Found", nropt, "solutions in", round(optimization_time, 2), "seconds\n")
    cat("Checked", points_checked, "/", length(search_order), "points (",
        round(100 * points_checked / length(search_order), 1), "% of search space)\n")
    cat("Best solution: rho =", round(rho_opt1, 3), ", psi =", round(psi_opt1, 3),
        ", ploidy =", round(ploidy_opt1, 3), ", goodness =", round(goodnessOfFit_opt1, 2), "\n")
  }

  # Generate plots (IDENTICAL to original)
  if (analysis=="paired"){
    if (!is.na(distancepng)) {
      png(filename = distancepng, width = 1000, height = 1000, res = 1000/7, type = "cairo")
    }
    ASCAT::ascat.plotSunrise(-d, psi_opt1_plot, rho_opt1_plot, minimise)
    if (!is.na(distancepng)) { dev.off() }
  }

  rho = rho_opt1
  psi = psi_opt1
  ploidy = ploidy_opt1

  nAfull = (rho-1-(b-1)*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
  nBfull = (rho-1+b*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
  nA = pmax(round(nAfull),0)
  nB = pmax(round(nBfull),0)

  rBacktransform = gamma*log((rho*(nA+nB)+(1-rho)*2)/((1-rho)*2+rho*psi),2)
  bBacktransform = (1-rho+rho*nB)/(2-2*rho+rho*(nA+nB))
  rConf = ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
  bConf = ifelse(bBacktransform!=0.5,pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))),NA)

  if(!is.na(reliabilityFile)){
    write.table(data.frame(segmentedBAF=b,backTransformedBAF=bBacktransform,confidenceBAF=bConf,segmentedR=r,backTransformedR=rBacktransform,confidenceR=rConf,nA=nA,nB=nB,nAfull=nAfull,nBfull=nBfull), reliabilityFile,sep=",",row.names=F)
  }
  confidence = ifelse(is.na(rConf),bConf,ifelse(is.na(bConf),rConf,(rConf+bConf)/2))

  # Create plots
  if (!is.na(copynumberprofilespng)) {
    png(filename = copynumberprofilespng, width = 2000, height = 500, res = 200, type = "cairo")
  }
  ASCAT::ascat.plotAscatProfile(n1all = nA, n2all = nB, heteroprobes = TRUE, ploidy = ploidy_opt1, rho = rho_opt1, goodnessOfFit = goodnessOfFit_opt1, nonaberrant = FALSE, ch = ch, lrr = lrr, bafsegmented = bafsegmented, chrs=chr.names)
  if (!is.na(copynumberprofilespng)) { dev.off() }

  if (!is.na(nonroundedprofilepng)) {
    png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200, type = "cairo")
  }
  ASCAT::ascat.plotNonRounded(ploidy = ploidy_opt1, rho = rho_opt1, goodnessOfFit = goodnessOfFit_opt1, nonaberrant = FALSE, nAfull = nAfull, nBfull = nBfull, bafsegmented = bafsegmented, ch = ch, lrr = lrr, chrs=chr.names)
  if (!is.na(nonroundedprofilepng)) { dev.off() }

  return(list(
    psi = psi,
    rho = rho,
    ploidy = ploidy,
    convergence_info = list(
      converged = TRUE,
      n_solutions_found = nropt,
      optimization_time = optimization_time,
      points_checked = points_checked,
      search_efficiency = points_checked / length(search_order)
    )
  ))
}

#' Create smart search order - best regions first
create_smart_search_order <- function(d, smart_ordering, verbose) {

  # Get all valid search points (excluding borders)
  nr <- nrow(d)
  nc <- ncol(d)
  search_points <- list()

  for (i in 4:(nr-3)) {
    for (j in 4:(nc-3)) {
      if (is.finite(d[i,j])) {
        search_points[[length(search_points) + 1]] <- list(i = i, j = j, distance = d[i,j])
      }
    }
  }

  if (!smart_ordering) {
    # Return in original order
    return(lapply(search_points, function(p) list(i = p$i, j = p$j)))
  }

  # Smart ordering: best distances first
  distances <- sapply(search_points, function(p) p$distance)
  order_idx <- order(distances)
  ordered_points <- search_points[order_idx]

  if (verbose) {
    cat("Smart ordering: searching best", length(ordered_points), "regions first\n")
    cat("Distance range:", round(min(distances), 4), "to", round(max(distances), 4), "\n")
  }

  return(lapply(ordered_points, function(p) list(i = p$i, j = p$j)))
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

#' Fast solution calculation (vectorized and optimized)
calculate_solution_fast <- function(psi, rho, s_b, s_r, s_length, total_length, gamma,
                                   min.ploidy, max.ploidy, min.rho, max.rho,
                                   min.goodness, distance_value, TheoretMaxdist, minimise,
                                   allow100percent, skip_zero_check = FALSE) {

  # Quick input validation
  if (is.na(psi) || is.na(rho) || psi <= 0 || rho <= 0 || rho > 1.1) {
    return(NULL)
  }

  # Quick constraint pre-check
  if (psi < min.ploidy || psi > max.ploidy || rho < min.rho || rho > max.rho) {
    return(NULL)
  }

  # Vectorized copy number calculation
  multiplier <- 2^(s_r/gamma) * ((1-rho)*2 + rho*psi)
  nA <- (rho - 1 - (s_b - 1) * multiplier) / rho
  nB <- (rho - 1 + s_b * multiplier) / rho

  # Quick validation
  if (any(is.na(nA)) || any(is.na(nB)) || any(!is.finite(nA)) || any(!is.finite(nB))) {
    return(NULL)
  }

  # Vectorized ploidy calculation
  ploidy <- sum((nA + nB) * s_length) / total_length

  if (is.na(ploidy) || !is.finite(ploidy) || ploidy <= 0) {
    return(NULL)
  }

  # Final ploidy constraint check
  if (ploidy < min.ploidy || ploidy > max.ploidy) {
    return(NULL)
  }

  # Fast goodness calculation
  if(minimise) {
    goodnessOfFit <- (1 - distance_value/TheoretMaxdist) * 100
  } else {
    goodnessOfFit <- -distance_value/TheoretMaxdist * 100
  }

  if (is.na(goodnessOfFit) || !is.finite(goodnessOfFit) || goodnessOfFit < min.goodness) {
    return(NULL)
  }

  # Zero check (only if needed)
  if (!skip_zero_check && !allow100percent) {
    nA_rounded <- round(nA)
    nB_rounded <- round(nB)

    percentzero <- (sum((nA_rounded == 0) * s_length) + sum((nB_rounded == 0) * s_length)) / total_length

    # Fast perczeroAbb calculation
    baf_mask <- s_b != 0.5
    if (any(baf_mask)) {
      denom <- sum(s_length[baf_mask])
      if (denom > 0) {
        perczeroAbb <- (sum((nA_rounded == 0) * s_length * baf_mask) +
                       sum((nB_rounded == 0) * s_length * baf_mask)) / denom
      } else {
        perczeroAbb <- 0
      }
    } else {
      perczeroAbb <- 0
    }

    if (is.na(perczeroAbb)) perczeroAbb <- 0

    if (!(percentzero > 0.01 || perczeroAbb > 0.1)) {
      return(NULL)
    }
  }

  return(list(
    psi = psi,
    rho = min(rho, 1.0),
    ploidy = ploidy,
    goodness = goodnessOfFit,
    distance = distance_value
  ))
}

#' Generate plots
generate_plots_battenberg <- function(analysis, distancepng, copynumberprofilespng, nonroundedprofilepng,
                                     d, psi_opt1, rho_opt1, ploidy_opt1, goodnessOfFit_opt1, minimise,
                                     b, r, s, gamma, ch, lrr, bafsegmented, chr.names, reliabilityFile) {
  
  if (analysis == "paired") {
    psi_opt1_plot <- psi_opt1
    rho_opt1_plot <- rho_opt1
    
    if (!is.na(distancepng)) {
      png(filename = distancepng, width = 1000, height = 1000, res = 1000/7, type = "cairo")
    }
    ASCAT::ascat.plotSunrise(-d, psi_opt1_plot, rho_opt1_plot, minimise)
    if (!is.na(distancepng)) { dev.off() }
  }
  
  nAfull <- (rho_opt1-1-(b-1)*2^(r/gamma)*((1-rho_opt1)*2+rho_opt1*psi_opt1))/rho_opt1
  nBfull <- (rho_opt1-1+b*2^(r/gamma)*((1-rho_opt1)*2+rho_opt1*psi_opt1))/rho_opt1
  nA <- pmax(round(nAfull),0)
  nB <- pmax(round(nBfull),0)
  
  if(!is.na(reliabilityFile)){
    rBacktransform <- gamma*log((rho_opt1*(nA+nB)+(1-rho_opt1)*2)/((1-rho_opt1)*2+rho_opt1*psi_opt1),2)
    bBacktransform <- (1-rho_opt1+rho_opt1*nB)/(2-2*rho_opt1+rho_opt1*(nA+nB))
    rConf <- ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
    bConf <- ifelse(bBacktransform!=0.5,pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))),NA)
    
    write.table(data.frame(segmentedBAF=b,backTransformedBAF=bBacktransform,confidenceBAF=bConf,
                          segmentedR=r,backTransformedR=rBacktransform,confidenceR=rConf,
                          nA=nA,nB=nB,nAfull=nAfull,nBfull=nBfull), 
                reliabilityFile, sep=",", row.names=F)
  }
  
  if (!is.na(copynumberprofilespng)) {
    png(filename = copynumberprofilespng, width = 2000, height = 500, res = 200, type = "cairo")
  }
  ASCAT::ascat.plotAscatProfile(n1all = nA, n2all = nB, heteroprobes = TRUE, 
                               ploidy = ploidy_opt1, rho = rho_opt1, goodnessOfFit = goodnessOfFit_opt1, 
                               nonaberrant = FALSE, ch = ch, lrr = lrr, bafsegmented = bafsegmented, 
                               chrs = chr.names)
  if (!is.na(copynumberprofilespng)) { dev.off() }
  
  if (!is.na(nonroundedprofilepng)) {
    png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200, type = "cairo")
  }
  ASCAT::ascat.plotNonRounded(ploidy = ploidy_opt1, rho = rho_opt1, goodnessOfFit = goodnessOfFit_opt1, 
                             nonaberrant = FALSE, nAfull = nAfull, nBfull = nBfull, 
                             bafsegmented = bafsegmented, ch = ch, lrr = lrr, chrs = chr.names)
  if (!is.na(nonroundedprofilepng)) { dev.off() }
}
