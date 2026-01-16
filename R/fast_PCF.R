# PCF-ALGORITHM (KL):
### EXACT version
exactPcf <- function(y, kmin = 5, gamma, yest) {
  ## Implementaion of exact PCF by Potts-filtering
  ## x: input array of (log2) copy numbers
  ## kmin: Mininal length of plateaus
  ## gamma: penalty for each discontinuity
  N <- length(y)
  yhat <- rep(0, N)
  if (N < 2 * kmin) {
    if (yest) {
      return(list(Lengde = N, sta = 1, mean = mean(y), nIntervals = 1, yhat = rep(mean(y), N)))
    } else {
      return(list(Lengde = N, sta = 1, mean = mean(y), nIntervals = 1))
    }
  }
  initSum <- sum(y[1:kmin])
  initKvad <- sum(y[1:kmin]^2)
  initAve <- initSum / kmin
  bestCost <- rep(0, N)
  bestCost[kmin] <- initKvad - initSum * initAve
  bestSplit <- rep(0, N)
  bestAver <- rep(0, N)
  bestAver[kmin] <- initAve
  Sum <- rep(0, N)
  Kvad <- rep(0, N)
  Aver <- rep(0, N)
  Cost <- rep(0, N)
  kminP1 <- kmin + 1
  for (k in (kminP1):(2 * kmin - 1)) {
    Sum[kminP1:k] <- Sum[kminP1:k] + y[k]
    Aver[kminP1:k] <- Sum[kminP1:k] / ((k - kmin):1)
    Kvad[kminP1:k] <- Kvad[kminP1:k] + y[k]^2
    bestAver[k] <- (initSum + Sum[kminP1]) / k
    bestCost[k] <- (initKvad + Kvad[kminP1]) - k * bestAver[k]^2
  }
  for (n in (2 * kmin):N) {
    yn <- y[n]
    yn2 <- yn^2
    Sum[kminP1:n] <- Sum[kminP1:n] + yn
    Aver[kminP1:n] <- Sum[kminP1:n] / ((n - kmin):1)
    Kvad[kminP1:n] <- Kvad[kminP1:n] + yn2
    nMkminP1 <- n - kmin + 1
    Cost[kminP1:nMkminP1] <- bestCost[kmin:(n - kmin)] + Kvad[kminP1:nMkminP1] - Sum[kminP1:nMkminP1] * Aver[kminP1:nMkminP1] + gamma
    Pos <- which.min(Cost[kminP1:nMkminP1]) + kmin
    cost <- Cost[Pos]
    aver <- Aver[Pos]
    totAver <- (Sum[kminP1] + initSum) / n
    totCost <- (Kvad[kminP1] + initKvad) - n * totAver * totAver
    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
      aver <- totAver
    }
    bestCost[n] <- cost
    bestAver[n] <- aver
    bestSplit[n] <- Pos - 1
  }
  n <- N
  antInt <- 0
  if (yest) {
    while (n > 0) {
      yhat[(bestSplit[n] + 1):n] <- bestAver[n]
      n <- bestSplit[n]
      antInt <- antInt + 1
    }
  } else {
    while (n > 0) {
      n <- bestSplit[n]
      antInt <- antInt + 1
    }
  }
  n <- N # nProbes   Spr Knut, fant ikke nProbes noe sted..
  lengde <- rep(0, antInt)
  start <- rep(0, antInt)
  verdi <- rep(0, antInt)
  oldSplit <- n
  antall <- antInt
  while (n > 0) {
    start[antall] <- bestSplit[n] + 1
    lengde[antall] <- oldSplit - bestSplit[n]
    verdi[antall] <- bestAver[n]
    n <- bestSplit[n]
    oldSplit <- n
    antall <- antall - 1
  }
  if (yest) {
    return(list(Lengde = lengde, sta = start, mean = verdi, nIntervals = antInt, yhat = yhat))
  } else {
    return(list(Lengde = lengde, sta = start, mean = verdi, nIntervals = antInt))
  }
}


selectFastPcf <- function(x, kmin, gamma, yest) {
  xLength <- length(x)
  if (xLength < 1000) {
    result <- runFastPcf(x, kmin, gamma, 0.15, 0.15, yest)
  } else {
    if (xLength < 15000) {
      result <- runFastPcf(x, kmin, gamma, 0.12, 0.05, yest)
    } else {
      result <- runPcfSubset(x, kmin, gamma, 0.12, 0.05, yest)
    }
  }
  return(result)
}


runFastPcf <- function(x, kmin, gamma, frac1, frac2, yest) {
  antGen <- length(x)
  mark <- filterMarkS4(x, kmin, 8, 1, frac1, frac2, 0.02, 0.9)
  mark[antGen] <- TRUE
  dense <- compact(x, mark)
  result <- PottsCompact(kmin, gamma, dense$Nr, dense$Sum, dense$Sq, yest)
  return(result)
}

runPcfSubset <- function(x, kmin, gamma, frac1, frac2, yest) {
  SUBSIZE <- 5000
  antGen <- length(x)
  mark <- filterMarkS4(x, kmin, 8, 1, frac1, frac2, 0.02, 0.9)
  markInit <- c(mark[1:(SUBSIZE - 1)], TRUE)
  compX <- compact(x[1:SUBSIZE], markInit)
  mark2 <- rep(FALSE, antGen)
  mark2[1:SUBSIZE] <- markWithPotts(kmin, gamma, compX$Nr, compX$Sum, compX$Sq, SUBSIZE)
  mark2[4 * SUBSIZE / 5] <- TRUE
  start <- 4 * SUBSIZE / 5 + 1
  while (start + SUBSIZE < antGen) {
    slutt <- start + SUBSIZE - 1
    markSub <- c(mark2[1:(start - 1)], mark[start:slutt])
    markSub[slutt] <- TRUE
    compX <- compact(x[1:slutt], markSub)
    mark2[1:slutt] <- markWithPotts(kmin, gamma, compX$Nr, compX$Sum, compX$Sq, slutt)
    start <- start + 4 * SUBSIZE / 5
    mark2[start - 1] <- TRUE
  }
  markSub <- c(mark2[1:(start - 1)], mark[start:antGen])
  compX <- compact(x, markSub)
  result <- PottsCompact(kmin, gamma, compX$Nr, compX$Sum, compX$Sq, yest)
  return(result)
}

PottsCompact <- function(kmin, gamma, nr, res, sq, yest) {
  ## Potts filtering on compact array;
  ## kmin: minimal length of plateau
  ## gamma: penalty for discontinuity
  ## nr: number of values between breakpoints
  ## res: sum of values between breakpoints
  ## sq: sum of squares of values between breakpoints

  N <- length(nr)
  Ant <- rep(0, N)
  Sum <- rep(0, N)
  Kvad <- rep(0, N)
  Cost <- rep(0, N)
  if (sum(nr) < 2 * kmin) {
    estim <- sum(res) / sum(nr)
    return(estim)
  }
  initAnt <- nr[1]
  initSum <- res[1]
  initKvad <- sq[1]
  initAve <- initSum / initAnt
  bestCost <- rep(0, N)
  bestCost[1] <- initKvad - initSum * initAve
  bestSplit <- rep(0, N)
  k <- 2
  while (sum(nr[1:k]) < 2 * kmin) {
    Ant[2:k] <- Ant[2:k] + nr[k]
    Sum[2:k] <- Sum[2:k] + res[k]
    Kvad[2:k] <- Kvad[2:k] + sq[k]
    bestCost[k] <- (initKvad + Kvad[2]) - (initSum + Sum[2])^2 / (initAnt + Ant[2])
    k <- k + 1
  }
  for (n in k:N) {
    Ant[2:n] <- Ant[2:n] + nr[n]
    Sum[2:n] <- Sum[2:n] + res[n]
    Kvad[2:n] <- Kvad[2:n] + sq[n]
    limit <- n
    while (limit > 2 && Ant[limit] < kmin) {
      limit <- limit - 1
    }
    Cost[2:limit] <- bestCost[1:limit - 1] + Kvad[2:limit] - Sum[2:limit]^2 / Ant[2:limit]
    Pos <- which.min(Cost[2:limit]) + 1
    cost <- Cost[Pos] + gamma
    totCost <- (Kvad[2] + initKvad) - (Sum[2] + initSum)^2 / (Ant[2] + initAnt)
    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
    }
    bestCost[n] <- cost
    bestSplit[n] <- Pos - 1
  }
  if (yest) {
    res <- findEst(bestSplit, N, nr, res, TRUE)
  } else {
    res <- findEst(bestSplit, N, nr, res, FALSE)
  }
  return(res)
}

compact <- function(y, mark) {
  ## accumulates numbers of observations, sums and
  ## sums of squares between potential breakpoints
  return(list(
    Nr = diff(append(0, which(mark))),
    Sum = diff(append(0, cumsum(y)[mark])),
    Sq = diff(append(0, cumsum(y^2)[mark]))
  ))
}

findEst <- function(bestSplit, N, Nr, Sum, yest) {
  n <- N
  lengde <- rep(0, N)
  antInt <- 0
  while (n > 0) {
    antInt <- antInt + 1
    lengde[antInt] <- n - bestSplit[n]
    n <- bestSplit[n]
  }
  lengde <- lengde[antInt:1]
  lengdeOrig <- rep(0, antInt)
  startOrig <- rep(1, antInt + 1)
  verdi <- rep(0, antInt)
  start <- rep(1, antInt + 1)
  for (i in 1:antInt) {
    start[i + 1] <- start[i] + lengde[i]
    lengdeOrig[i] <- sum(Nr[start[i]:(start[i + 1] - 1)])
    startOrig[i + 1] <- startOrig[i] + lengdeOrig[i]
    verdi[i] <- sum(Sum[start[i]:(start[i + 1] - 1)]) / lengdeOrig[i]
  }

  if (yest) {
    yhat <- rep(0, startOrig[antInt + 1] - 1)
    for (i in 1:antInt) {
      yhat[startOrig[i]:(startOrig[i + 1] - 1)] <- verdi[i]
    }
    startOrig <- startOrig[1:antInt]
    return(list(Lengde = lengdeOrig, sta = startOrig, mean = verdi, nIntervals = antInt, yhat = yhat))
  } else {
    startOrig <- startOrig[1:antInt]
    return(list(Lengde = lengdeOrig, sta = startOrig, mean = verdi, nIntervals = antInt))
  }
}


markWithPotts <- function(kmin, gamma, nr, res, sq, subsize) {
  ## Potts filtering on compact array;
  ## kmin: minimal length of plateau
  ## gamma: penalty for discontinuity
  ## nr: number of values between breakpoints
  ## res: sum of values between breakpoints
  ## sq: sum of squares of values between breakpoints

  N <- length(nr)
  Ant <- rep(0, N)
  Sum <- rep(0, N)
  Kvad <- rep(0, N)
  Cost <- rep(0, N)
  markSub <- rep(FALSE, N)
  initAnt <- nr[1]
  initSum <- res[1]
  initKvad <- sq[1]
  initAve <- initSum / initAnt
  bestCost <- rep(0, N)
  bestCost[1] <- initKvad - initSum * initAve
  bestSplit <- rep(0, N)
  k <- 2
  while (sum(nr[1:k]) < 2 * kmin) {
    Ant[2:k] <- Ant[2:k] + nr[k]
    Sum[2:k] <- Sum[2:k] + res[k]
    Kvad[2:k] <- Kvad[2:k] + sq[k]
    bestCost[k] <- (initKvad + Kvad[2]) - (initSum + Sum[2])^2 / (initAnt + Ant[2])
    k <- k + 1
  }
  for (n in k:N) {
    Ant[2:n] <- Ant[2:n] + nr[n]
    Sum[2:n] <- Sum[2:n] + res[n]
    Kvad[2:n] <- Kvad[2:n] + sq[n]
    limit <- n
    while (limit > 2 && Ant[limit] < kmin) {
      limit <- limit - 1
    }
    Cost[2:limit] <- bestCost[1:limit - 1] + Kvad[2:limit] - Sum[2:limit]^2 / Ant[2:limit]
    Pos <- which.min(Cost[2:limit]) + 1
    cost <- Cost[Pos] + gamma
    totCost <- (Kvad[2] + initKvad) - (Sum[2] + initSum)^2 / (Ant[2] + initAnt)
    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
    }
    bestCost[n] <- cost
    bestSplit[n] <- Pos - 1
    markSub[Pos - 1] <- TRUE
  }
  help <- findMarks(markSub, nr, subsize)
  return(help = help)
}


findMarks <- function(markSub, Nr, subsize) {
  ## markSub: marks in compressed scale
  ## NR: number of observations between potenstial breakpoints
  mark <- rep(FALSE, subsize) ## marks in original scale
  if (sum(markSub) < 1) {
    return(mark)
  } else {
    N <- length(markSub)
    ant <- seq(1:N)
    help <- ant[markSub]
    lengdeHelp <- length(help)
    help0 <- c(0, help[1:(lengdeHelp - 1)])
    lengde <- help - help0
    start <- 1
    oldStart <- 1
    startOrig <- 1
    for (i in 1:lengdeHelp) {
      start <- start + lengde[i]
      lengdeOrig <- sum(Nr[oldStart:(start - 1)])
      startOrig <- startOrig + lengdeOrig
      mark[startOrig - 1] <- TRUE
      oldStart <- start
    }
    return(mark)
  }
}

filterMarkS4 <- function(x, kmin, L, L2, frac1, frac2, frac3, thres) {
  lengdeArr <- length(x)
  xc <- c(0, cumsum(x)) # Lead with 0 so xc[1] is 0

  # --- Cost 1 Calculation (Window L) ---
  ind11 <- 1:(lengdeArr - 6 * L + 1)
  # The formula: 4*xc[ind13] - xc[ind11] - xc[ind12] - xc[ind14] - xc[ind15]
  cost1 <- abs(4 * xc[ind11 + 3 * L] - xc[ind11] - xc[ind11 + L] - xc[ind11 + 5 * L] - xc[ind11 + 6 * L])
  cost1_full <- c(numeric(3 * L - 1), cost1, numeric(3 * L))

  # --- Rolling Max Parity ---
  # Your pmax was: pmax(cost1[i], cost1[i+1], ..., cost1[i+6])
  # To match 'rep(0, 3)' at both ends, we use align="center" with a window of 7
  test1 <- RcppRoll::roll_max(cost1_full, n = 7, fill = 0, align = "center")

  cost1B <- cost1_full[cost1_full >= thres * test1]
  frac1B <- min(0.8, frac1 * length(cost1_full) / length(cost1B))
  limit1 <- collapse::fquantile(cost1B, (1 - frac1B), names = FALSE)
  mark <- (cost1_full > limit1) & (cost1_full > 0.9 * test1)

  # --- Cost 2 Calculation (Window L2) ---
  ind21 <- 1:(lengdeArr - 6 * L2 + 1)
  cost2 <- abs(4 * xc[ind21 + 3 * L2] - xc[ind21] - xc[ind21 + L2] - xc[ind21 + 5 * L2] - xc[ind21 + 6 * L2])
  limit2 <- collapse::fquantile(cost2, (1 - frac2), names = FALSE)

  mark2_core <- (cost2 > limit2)
  mark2 <- c(numeric(3 * L2 - 1), mark2_core, numeric(3 * L2))

  # --- Edge Case Overrides ---
  if (3 * L > kmin) {
    mark[kmin:(3 * L - 1)] <- TRUE
    mark[(lengdeArr - 3 * L + 1):(lengdeArr - kmin)] <- TRUE
  } else {
    mark[kmin] <- TRUE
    mark[lengdeArr - kmin] <- TRUE
  }

  # --- Short Segment Detection (kmin) ---
  if (kmin > 1) {
    i_s <- 1:(lengdeArr - 3 * kmin + 1)
    shortAb <- abs(3 * (xc[i_s + 2 * kmin] - xc[i_s + kmin]) - (xc[i_s + 3 * kmin] - xc[i_s]))

    test_s <- RcppRoll::roll_max(shortAb, n = 7, fill = 0, align = "center")

    cost1C <- shortAb[shortAb >= thres * test_s]
    frac1C <- min(0.8, frac3 * length(shortAb) / length(cost1C))
    limit3 <- collapse::fquantile(cost1C, (1 - frac1C), names = FALSE)

    markH1 <- (shortAb > limit3) & (shortAb > thres * test_s)

    # Pixel-perfect shift reproduction
    markH2 <- c(logical(kmin - 1), markH1, logical(2 * kmin))
    markH3 <- c(logical(2 * kmin - 1), markH1, logical(kmin))
    mark <- mark | mark2 | markH2 | markH3
  } else {
    mark <- mark | mark2
  }

  # --- Final Boundary Cleanup ---
  # Re-applying the final mark overrides exactly as the original function
  if (3 * L > kmin) {
    mark[1:(kmin - 1)] <- FALSE
    mark[kmin:(3 * L - 1)] <- TRUE
    mark[(lengdeArr - 3 * L + 1):(lengdeArr - kmin)] <- TRUE
    mark[(lengdeArr - kmin + 1):(lengdeArr - 1)] <- FALSE
  } else {
    mark[1:(kmin - 1)] <- FALSE
    mark[(lengdeArr - kmin + 1):(lengdeArr - 1)] <- FALSE
    mark[kmin] <- TRUE
    mark[lengdeArr - kmin] <- TRUE
  }
  mark[lengdeArr] <- TRUE

  return(mark)
}

# Optimized function to calculate the Median Absolute Deviation of a signal
# after removing a running median trend.
get_mad <- function(x, k = 25) {
  # Use collapse for fast, memory-efficient subsetting
  # Removes zeros which often represent missing/imputed data in genomics
  x_filtered <- collapse::fsubset(x, x != 0)

  # Use rlang to safely check for empty input after filtering
  if (length(x_filtered) == 0) {
    return(NA)
  }

  # Calculate running median parameters
  n <- length(x_filtered)
  filt_width <- 2 * k + 1

  # Ensure filt_width is odd and does not exceed n to satisfy runmed requirements
  if (filt_width > n) {
    filt_width <- if (n %% 2 == 0) max(1, n - 1) else max(1, n)
  }

  # Calculate the running median using the C-based engine
  # endrule = "median" ensures we don't get NAs at the start/end of the vector
  run_median <- stats::runmed(x_filtered, k = filt_width, endrule = "median")

  # Calculate the difference and the MAD
  # collapse::fmad is significantly faster than stats::mad
  residual_signal <- x_filtered - run_median
  SD <- stats::mad(residual_signal)

  return(SD)
}
