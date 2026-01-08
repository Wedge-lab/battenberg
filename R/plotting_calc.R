########################################################################################
# Various functions for calculating from data for plotting
########################################################################################
#' Calc copy number of major allele per segment from a subclones data.frame
#' @noRd
calc_total_cn_major <- function(bb) {
  return(bb$nMaj1_A * bb$frac1_A + ifelse(bb$frac1_A < 1, bb$nMaj2_A * bb$frac2_A, 0))
}

#' Calc copy number of minor allele per segment from a subclones data.frame
#' @noRd
calc_total_cn_minor <- function(bb) {
  return(bb$nMin1_A * bb$frac1_A + ifelse(bb$frac1_A < 1, bb$nMin2_A * bb$frac2_A, 0))
}

#' Calc total copy number per segment from a subclones data.frame
#' @noRd
calculate_bb_total_cn <- function(bb) {
  return((bb$nMaj1_A + bb$nMin1_A) * bb$frac1_A + ifelse(!is.na(bb$frac2_A), (bb$nMaj2_A + bb$nMin2_A) * bb$frac2_A, 0))
}

#' Calc ploidy from a subclones data.frame
#' @noRd
calc_ploidy <- function(bb) {
  bb$len <- bb$endpos / 1000 - bb$startpos / 1000
  bb$total_cn <- calculate_bb_total_cn(bb)
  ploidy <- sum(bb$total_cn * bb$len) / sum(bb$len)
  return(ploidy)
}

#' Transform logR into an estimate of total copy number given purity and total ploidy (tumour+normal)
#' @noRd
logr2tumcn <- function(cellularity, total_ploidy, logR) {
  return(((total_ploidy * (2^logR)) - 2 * (1 - cellularity)) / cellularity)
}

#' Calc psi from psi_t and rho
#' @noRd
psit2psi <- function(rho, psi_t) {
  return(rho * psi_t + 2 * (1 - rho))
}

#' Calc psi_t from psi and rho
#' @noRd
psi2psit <- function(rho, psi) {
  return((psi - 2 * (1 - rho)) / rho)
}
