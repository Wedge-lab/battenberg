#' Generate plots
generate_plots_battenberg <- function(
  analysis, distancepng, copynumberprofilespng, nonroundedprofilepng,
  d, psi_opt1, rho_opt1, ploidy_opt1, goodness_of_fit_opt1, minimise,
  b, r, s, gamma, ch, lrr, bafsegmented, chr_names, reliabilityFile
) {
  if (analysis == "paired") {
    psi_opt1_plot <- psi_opt1
    rho_opt1_plot <- rho_opt1

    if (!is.na(distancepng)) {
      grDevices::png(filename = distancepng, width = 1000, height = 1000, res = 1000 / 7, type = "cairo")
    }
    ASCAT::ascat.plotSunrise(-d, psi_opt1_plot, rho_opt1_plot, minimise)
    if (!is.na(distancepng)) {
      grDevices::dev.off()
    }
  }

  nAfull <- (rho_opt1 - 1 - (b - 1) * 2^(r / gamma) * ((1 - rho_opt1) * 2 + rho_opt1 * psi_opt1)) / rho_opt1
  nBfull <- (rho_opt1 - 1 + b * 2^(r / gamma) * ((1 - rho_opt1) * 2 + rho_opt1 * psi_opt1)) / rho_opt1
  nA <- pmax(round(nAfull), 0)
  nB <- pmax(round(nBfull), 0)

  if (!is.na(reliabilityFile)) {
    rBacktransform <- gamma * log((rho_opt1 * (nA + nB) + (1 - rho_opt1) * 2) / ((1 - rho_opt1) * 2 + rho_opt1 * psi_opt1), 2)
    bBacktransform <- (1 - rho_opt1 + rho_opt1 * nB) / (2 - 2 * rho_opt1 + rho_opt1 * (nA + nB))
    rConf <- ifelse(abs(rBacktransform) > 0.15, pmin(100, pmax(0, 100 * (1 - abs(rBacktransform - r) / abs(r)))), NA)
    bConf <- ifelse(bBacktransform != 0.5, pmin(100, pmax(0, ifelse(b == 0.5, 100, 100 * (1 - abs(bBacktransform - b) / abs(b - 0.5))))), NA)

    data.table::fwrite(
      data.frame(
        segmentedBAF = b, backTransformedBAF = bBacktransform, confidenceBAF = bConf,
        segmentedR = r, backTransformedR = rBacktransform, confidenceR = rConf,
        nA = nA, nB = nB, nAfull = nAfull, nBfull = nBfull
      ),
      reliabilityFile,
      sep = ",", row.names = F
    )
  }

  if (!is.na(copynumberprofilespng)) {
    grDevices::png(
      filename = copynumberprofilespng,
      width = 2000, height = 500,
      res = 200, type = "cairo"
    )
  }
  ASCAT::ascat.plotAscatProfile(
    n1all = nA, n2all = nB, heteroprobes = TRUE,
    ploidy = ploidy_opt1, rho = rho_opt1,
    goodness_of_fit = goodness_of_fit_opt1,
    nonaberrant = FALSE, ch = ch, lrr = lrr, bafsegmented = bafsegmented,
    chrs = chr_names
  )
  if (!is.na(copynumberprofilespng)) {
    grDevices::dev.off()
  }

  if (!is.na(nonroundedprofilepng)) {
    grDevices::png(
      filename = nonroundedprofilepng,
      width = 2000, height = 500,
      res = 200, type = "cairo"
    )
  }
  ASCAT::ascat.plotNonRounded(
    ploidy = ploidy_opt1, rho = rho_opt1,
    goodness_of_fit = goodness_of_fit_opt1,
    nonaberrant = FALSE, nAfull = nAfull,
    nBfull = nBfull,
    bafsegmented = bafsegmented, ch = ch, lrr = lrr, chrs = chr_names
  )
  if (!is.na(nonroundedprofilepng)) {
    grDevices::dev.off()
  }
}
