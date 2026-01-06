.onLoad <- function(libname, pkgname) {
  # Keep your scipen setting
  logger::log_threshold(logger::INFO, namespace = pkgname)
  options(scipen = 999)
}
