#' @importFrom gtools mixedsort
NULL

#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
create_haplotype_plot <- function(
  chrom_position,
  points.blue, points.red,
  x_min, x_max,
  title, xlab, ylab
) {
  graphics::par(
    pch = ".", cex = 1, cex.main = 0.8,
    cex.axis = 0.6, cex.lab = 0.7,
    yaxp = c(-0.05, 1.05, 6)
  )
  graphics::plot(
    c(x_min, x_max), c(0, 1),
    type = "n", main = title, xlab = xlab, ylab = ylab
  )
  if (length(chrom_position) > 0) {
    graphics::points(chrom_position, points.blue, col = "blue")
    graphics::points(chrom_position, points.red, col = "red")
  }
}

#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
create_segmented_plot <- function(
  chrom_position,
  points.red,
  points.green,
  x_min, x_max,
  title, xlab,
  ylab,
  prior_bkps_pos = NULL
) {
  graphics::par(
    mar = c(5, 5, 5, 0.5),
    cex = 0.4,
    cex.main = 3,
    cex.axis = 2,
    cex.lab = 2
  )
  graphics::plot(
    c(x_min, x_max), c(0, 1),
    pch = ".", type = "n", main = title, xlab = xlab, ylab = ylab
  )
  graphics::points(
    chrom_position, points.red,
    pch = ".", col = "red", cex = 2
  )
  graphics::points(
    chrom_position, points.green,
    pch = 19, cex = 0.5, col = "green"
  )
  if (!is.null(prior_bkps_pos)) {
    for (i in seq_along(prior_bkps_pos)) {
      graphics::abline(v = prior_bkps_pos[i])
    }
  }
}

#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
create_baf_plot <- function(
  chrom_position,
  points_red_blue, plot_red,
  points_darkred, points_darkblue,
  x_min, x_max,
  title, xlab, ylab,
  prior_bkps_pos = NULL
) {
  graphics::par(
    mar = c(5, 5, 5, 0.5), cex = 0.4, cex.main = 3, cex.axis = 2, cex.lab = 2
  )
  graphics::plot(
    c(x_min, x_max), c(0, 1),
    pch = ".", type = "n", main = title, xlab = xlab, ylab = ylab
  )
  graphics::points(
    chrom_position, points_red_blue,
    pch = ".", col = ifelse(plot_red, "red", "blue"), cex = 2
  )
  graphics::points(
    chrom_position, points_darkred,
    pch = 19, cex = 0.5, col = "darkred"
  )
  graphics::points(
    chrom_position, points_darkblue,
    pch = 19, cex = 0.5, col = "darkblue"
  )
  if (!is.null(prior_bkps_pos)) {
    for (i in seq_along(prior_bkps_pos)) {
      graphics::abline(v = prior_bkps_pos[i])
    }
  }
}

#' Function that creates the plots for subclonal copy number
#' Note: This is a plot PER chromosome.
#' @noRd
create_subclonal_cn_plot <- function(
  chrom,
  chrom_position,
  LogRposke, LogRchr,
  BAFchr, BAFsegchr,
  BAFpvalschr, subcloneres,
  siglevel, x_min, x_max,
  title, xlab, ylab_logr,
  ylab_baf, breakpoints_pos = NULL,
  svs_pos = NULL
) {
  plot_breakpoints <- function(breakpoints, svs_pos) {
    # Plot the breakpoints
    if (!is.null(breakpoints)) {
      for (i in seq_along(breakpoints)) {
        graphics::abline(v = breakpoints[i], col = "darkgrey", lwd = 1)
      }
    }

    # Overplot the SV breakpoints, if supplied
    if (!is.null(svs_pos)) {
      for (i in seq_along(svs_pos)) {
        graphics::abline(v = svs_pos[i], lty = 3, col = "lightgreen", lwd = 1)
      }
    }
  }

  # Plot the logR
  graphics::par(
    mar = c(2.5, 2.5, 2.5, 0.25),
    cex = 0.4, cex.main = 1.5,
    cex.axis = 1, cex.lab = 1, mfrow = c(2, 1)
  )
  graphics::plot(
    c(x_min, x_max), c(-3, 3),
    pch = ".", type = "n",
    main = title, xlab = xlab,
    ylab = ylab_logr
  )
  graphics::points(LogRposke / 1000000, LogRchr, pch = ".", col = "grey")
  plot_breakpoints(breakpoints_pos, svs_pos)

  # Plot BAF
  graphics::plot(
    c(x_min, x_max),
    c(0, 1),
    pch = ".", type = "n",
    main = title,
    xlab = xlab, ylab = ylab_baf
  )
  graphics::points(chrom_position, BAFchr, pch = ".", col = "grey")
  plot_breakpoints(breakpoints_pos, svs_pos)

  # Plot segments in top of BAF
  graphics::points(
    chrom_position, BAFsegchr,
    pch = 19, cex = 0.5,
    col = ifelse(BAFpvalschr > siglevel,
      "darkgreen", "red"
    )
  )
  graphics::points(
    chrom_position, 1 - BAFsegchr,
    pch = 19, cex = 0.5,
    col = ifelse(BAFpvalschr > siglevel,
      "darkgreen", "red"
    )
  )
  for (i in seq_len(dim(subcloneres)[1])) {
    if (subcloneres[i, 1] == chrom) {
      graphics::text(
        (as.numeric(subcloneres[i, "startpos"]) + as.numeric(subcloneres[i, "endpos"])) / 2 / 1000000, as.numeric(subcloneres[i, "BAF"]) - 0.04,
        paste(subcloneres[i, "nMaj1_A"], "+", subcloneres[i, "nMin1_A"], ": ", 100 * round(as.numeric(subcloneres[i, "frac1_A"]), 3), "%", sep = ""),
        cex = 0.8
      )
      if (!is.na(subcloneres[i, "nMaj2_A"])) {
        graphics::text((as.numeric(subcloneres[i, "startpos"]) + as.numeric(subcloneres[i, "endpos"])) / 2 / 1000000, as.numeric(subcloneres[i, "BAF"]) - 0.08,
          paste(subcloneres[i, "nMaj2_A"], "+", subcloneres[i, "nMin2_A"], ": ", 100 * round(as.numeric(subcloneres[i, "frac2_A"]), 3), "%", sep = ""),
          cex = 0.8
        )
      }
    }
  }
}


#' Make GW CN plot where subclonal CN is represented as a mixture of two states
#' NAP - July 2020 - updated main title now replacing 'cellularity' with 'purity' and 'goodness-of-fit' with 'PGAclonal' + adding TUMOURNAME
#' NAP - November 2023 - Replacing 'PGAclonal' with 'PGA.is.clonal' for more clarity
#' @noRd
create_bb_plot_average <- function(
  bafsegmented, ploidy, rho,
  goodness_of_fit, pos_min, pos_max,
  segment_states_min, segment_states_tot,
  chr_segs, chr_names, tumourname, ylim = 5
) {
  # Plot main frame and title
  graphics::par(
    mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main = 3, cex.axis = 2.5
  )
  maintitle <- paste0(
    substring(
      tumourname, 36,
      first = TRUE
    ),
    ", Ploidy: ", sprintf("%1.2f", ploidy),
    ", Purity: ", sprintf("%2.0f", rho * 100),
    "%, PGA.is.clonal: ",
    sprintf("%2.1f", goodness_of_fit * 100), "%"
  )
  graphics::plot(
    c(1, nrow(bafsegmented)), c(0, ylim),
    type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = ""
  )
  graphics::abline(v = 0, lty = 1, col = "lightgrey")
  # Horizontal lines for y=0 to y=5
  graphics::abline(h = c(0:ylim), lty = 1, col = "lightgrey")
  # Minor allele in gray, total CN in orange
  graphics::segments(
    x0 = pos_min,
    y0 = segment_states_min,
    x1 = pos_max,
    y1 = segment_states_min,
    col = "#2f4f4f", pch = "|", lwd = 6, lend = 1
  )
  graphics::segments(
    x0 = pos_min, y0 = segment_states_tot,
    x1 = pos_max, y1 = segment_states_tot,
    col = "#E69F00", pch = "|", lwd = 6, lend = 1
  )

  # Plot the vertical lines that show start/end of a chromosome
  chrk_tot_len <- 0
  for (i in seq_along(chr_segs)) {
    chrk <- chr_segs[[i]]
    chrk_hetero <- names(bafsegmented)[chrk]
    chrk_tot_len_prev <- chrk_tot_len
    chrk_tot_len <- chrk_tot_len + length(chrk_hetero)
    vpos <- chrk_tot_len
    tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
    graphics::text(tpos, ylim, chr_names[i], pos = 1, cex = 2)
    graphics::abline(v = vpos, lty = 1, col = "lightgrey")
  }
}

#' Make GW CN plot where subclonal CN is represented by two separate states
#' NAP - July 2020 - updated main title now replacing 'cellularity' with 'purity' and 'goodness-of-fit' with 'PGAclonal' + adding TUMOURNAME
#' NAP - November 2023 - Replacing 'PGAclonal' with 'PGA.is.clonal' for more clarity
#' @noRd
create_bb_plot_subclones <- function(
  bafsegmented, subclones, ploidy,
  rho, goodness_of_fit, pos_min,
  pos_max, subcl_min, subcl_max,
  is_subclonal, is_subclonal_maj,
  is_subclonal_min, chr_segs,
  chr_names, tumourname, ylim = 5
) {
  graphics::par(
    mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main = 3, cex.axis = 2.5
  )
  maintitle <- paste0(
    substring(tumourname, 36, first = TRUE),
    ", Ploidy: ", sprintf("%1.2f", ploidy),
    ", Purity: ", sprintf("%2.0f", rho * 100),
    "%, PGA.is.clonal: ",
    sprintf("%2.1f", goodness_of_fit * 100), "%"
  )

  graphics::plot(
    c(1, nrow(bafsegmented)), c(0, ylim),
    type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = ""
  )
  graphics::abline(
    v = 0, lty = 1, col = "lightgrey"
  )
  # Minor allele clonal and lowest of the two states when subclonal
  graphics::segments(
    x0 = pos_min, y0 = subclones$nMin1_A - 0.1,
    x1 = pos_max, y1 = subclones$nMin1_A - 0.1, col = "#2f4f4f", pch = "|",
    lwd = ifelse(is_subclonal_min, 6 * subclones$frac1_A, 6), lend = 1
  )

  if (sum(is_subclonal) > 0) {
    # Minor allele highest of the two states when subclonal
    graphics::segments(
      x0 = subcl_min, y0 = subclones$nMin2_A[is_subclonal] - 0.1,
      x1 = subcl_max, y1 = subclones$nMin2_A[is_subclonal] - 0.1, col = "#2f4f4f", pch = "|",
      lwd = ifelse(
        is_subclonal_min[is_subclonal],
        6 * subclones$frac2_A[is_subclonal], 0
      ), lend = 1
    )

    # Total CN, when minor allele subclonal CN (one of the two alleles)
    graphics::segments(
      x0 = subcl_min, y0 = subclones$nMaj1_A[is_subclonal] + subclones$nMin1_A[is_subclonal] + 0.1,
      x1 = subcl_max, y1 = subclones$nMaj1_A[is_subclonal] + subclones$nMin1_A[is_subclonal] + 0.1, col = "#E69F00", pch = "|",
      lwd = ifelse(
        is_subclonal_min[is_subclonal], 6 * subclones$frac1_A[is_subclonal], 0
      ), lend = 1
    )

    # Total CN, when minor allele subclonal CN (the other allele)
    graphics::segments(
      x0 = subcl_min, y0 = subclones$nMaj2_A[is_subclonal] + subclones$nMin2_A[is_subclonal] + 0.1,
      x1 = subcl_max, y1 = subclones$nMaj2_A[is_subclonal] + subclones$nMin2_A[is_subclonal] + 0.1, col = "#E69F00", pch = "|",
      lwd = ifelse(
        is_subclonal_min[is_subclonal], 6 * subclones$frac2_A[is_subclonal], 0
      ), lend = 1
    )
  }

  # Total CN, when major allele clonal and subclonal, unless the minor allele is subclonal (then plot nothing, done above)
  graphics::segments(
    x0 = pos_min, y0 = subclones$nMaj1_A + subclones$nMin1_A + 0.1,
    x1 = pos_max, y1 = subclones$nMaj1_A + subclones$nMin1_A + 0.1, col = "#E69F00", pch = "|",
    lwd = ifelse(
      is_subclonal_maj & (!is_subclonal_min), 6 * subclones$frac1_A, 0
    ), lend = 1
  )

  # Total CN, when subclonal major allele and not subclonal minor allele (the other allele)
  graphics::segments(
    x0 = pos_min, y0 = subclones$nMaj2_A + subclones$nMin2_A + 0.1,
    x1 = pos_max, y1 = subclones$nMaj2_A + subclones$nMin2_A + 0.1, col = "#E69F00", pch = "|",
    lwd = ifelse(
      is_subclonal_maj & (!is_subclonal_min), 6 * subclones$frac2_A, 0
    ), lend = 1
  )

  # Total allele when major and minor both non-subclonal
  graphics::segments(
    x0 = pos_min, y0 = subclones$nMaj1_A + subclones$nMin1_A + 0.1,
    x1 = pos_max, y1 = subclones$nMaj1_A + subclones$nMin1_A + 0.1, col = "#E69F00", pch = "|",
    lwd = ifelse((!is_subclonal_maj) & (!is_subclonal_min), 6, 0), lend = 1
  )

  chrk_tot_len <- 0
  for (i in seq_along(chr_segs)) {
    chrk <- chr_segs[[i]]
    chrk_hetero <- names(bafsegmented)[chrk]
    chrk_tot_len_prev <- chrk_tot_len
    chrk_tot_len <- chrk_tot_len + length(chrk_hetero)
    vpos <- chrk_tot_len
    tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
    graphics::text(tpos, ylim, chr_names[i], pos = 1, cex = 2)
    graphics::abline(v = vpos, lty = 1, col = "lightgrey")
  }
}

#' Code extracted from the plot in clonal_ascat find_centroid_of_global_minima.
#' Note: This is a temporary function and VERY similar to clonal_runascat.plot1()
#' @noRd
clonal_findcentroid_plot <- function(minimise, dist_choice, d, psis, rhos, new_bounds) {
  graphics::par(
    mar = c(5, 5, 0.5, 0.5), cex = 0.75, cex.lab = 2, cex.axis = 2
  )
  # DCW 240314 reverse colour palette, so blue always corresponds to best region
  if (minimise) {
    hmcol <- rev(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(10, "RdBu")
      )(256)
    )
  } else {
    hmcol <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(10, "RdBu")
    )(256)
  }
  if (dist_choice == 4) {
    graphics::image(
      d,
      col = hmcol, axes = FALSE,
      xlab = "Ploidy", ylab = "Aberrant cell fraction"
    )
  } else {
    graphics::image(
      log(d),
      col = hmcol, axes = FALSE,
      xlab = "Ploidy", ylab = "Aberrant cell fraction"
    )
  }
  psi_min <- new_bounds$psi_min
  psi_max <- new_bounds$psi_max
  rho_min <- new_bounds$rho_min
  rho_max <- new_bounds$rho_max

  psi_range <- psi_max - psi_min
  rho_range <- rho_max - rho_min

  psi_min_label <- ceiling(10 * psi_min) / 10
  psi_max_label <- floor(10 * psi_max) / 10
  psi_label_interval <- 0.1

  psi_min_label_standardised <- (psi_min_label - psi_min) / psi_range
  psi_max_label_standardised <- (psi_max_label - psi_min) / psi_range
  psi_label_interval_standardised <- psi_label_interval / psi_range

  rho_min_label <- ceiling(100 * rho_min) / 100
  rho_max_label <- floor(100 * rho_max) / 100
  rho_label_interval <- 0.01

  rho_min_label_standardised <- (rho_min_label - rho_min) / rho_range
  rho_max_label_standardised <- (rho_max_label - rho_min) / rho_range
  rho_label_interval_standardised <- rho_label_interval / rho_range

  graphics::axis(
    1,
    at = seq(
      psi_min_label_standardised,
      psi_max_label_standardised,
      by = psi_label_interval_standardised
    ),
    labels = seq(psi_min_label, psi_max_label, by = psi_label_interval)
  )
  graphics::axis(
    2,
    at = seq(
      rho_min_label_standardised,
      rho_max_label_standardised,
      by = rho_label_interval_standardised
    ),
    labels = seq(rho_min_label, rho_max_label, by = rho_label_interval)
  )

  graphics::points(
    (psis - psi_min) / psi_range, (rhos - rho_min) / rho_range,
    col = c("green", "darkgreen"), pch = "X", cex = 2
  )
}

# Plot Battenberg copy number solutions for a segment
# Refactored for clarity and data.table integration
squaresplot <- function(tumourname, run_dir, segment_chr, segment_pos,
                        platform_gamma = 1, pdf = 0, binwidth_baf = 0.25, xylimits = c(-0.2, 5)) {
  # Construct output paths
  ext <- if (pdf) ".pdf" else ".png"
  out_file <- file.path(run_dir, paste0(tumourname, "_squares_chr", segment_chr, "_", segment_pos, ext))

  if (pdf) {
    grDevices::pdf(file = out_file, width = 7, height = 7)
  } else {
    grDevices::png(filename = out_file, width = 1200, height = 1200, res = 200, type = "cairo")
  }

  # Parse chromosomal position
  segment_pos_num <- as.numeric(gsub("M", "000000", segment_pos))

  # Read data using data.table
  cn_file <- file.path(run_dir, paste0(tumourname, "_copynumber.txt"))
  subclones <- data.table::fread(cn_file, data.table = FALSE)

  # Select specific segment
  subclone <- subclones[(subclones$chr == segment_chr) &
    (subclones$startpos <= segment_pos_num) &
    (subclones$endpos >= segment_pos_num), ]

  # Get best rho and psi parameters
  rp_file <- file.path(run_dir, paste0(tumourname, "_rho_and_psi.txt"))
  rhopsi_df <- data.table::fread(rp_file, data.table = FALSE)
  rhopsi <- rhopsi_df[rhopsi_df$is_best == TRUE, c("rho", "psi")]

  rho <- rhopsi$rho
  psi <- rhopsi$psi

  # Theoretical calculations
  logr_comp <- 2^(subclone$LogR / platform_gamma)
  p_comp <- ((1 - rho) * 2 + rho * psi)
  nMincalc <- (rho - 1 - (subclone$BAF - 1) * logr_comp * p_comp) / rho
  nMajcalc <- (rho - 1 + subclone$BAF * logr_comp * p_comp) / rho

  # Grid function
  isobafline <- function(nB, cstbaf) {
    (1 - rho + rho * nB - cstbaf * (2 - 2 * rho) - rho * cstbaf * nB) / (rho * cstbaf)
  }

  # Base Plot
  q <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(name = "nMajor", breaks = 0:max(xylimits), limits = xylimits) +
    ggplot2::scale_y_continuous(name = "nMinor", breaks = 0:max(xylimits), limits = xylimits) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "darkgrey", size = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Grid Lines
  baf_seq <- seq(0, 1, binwidth_baf)
  for (bafval in baf_seq) {
    q <- q + ggplot2::stat_function(fun = isobafline, args = list(cstbaf = bafval), colour = "blue", alpha = 0.6)
  }
  q <- q + ggplot2::stat_function(fun = isobafline, args = list(cstbaf = subclone$BAF), colour = "green")

  # Isologrline (Red Segment)
  err_df <- data.frame(
    x = floor(nMajcalc) - 0.2, y = ceiling(nMincalc) + 0.2,
    xend = ceiling(nMajcalc) + 0.2, yend = floor(nMincalc) - 0.2
  )

  # Note the use of rlang::.data here
  q <- q + ggplot2::geom_segment(
    data = err_df,
    ggplot2::aes(
      x = rlang::.data$x,
      y = rlang::.data$y,
      xend = rlang::.data$xend,
      yend = rlang::.data$yend
    ),
    colour = "red", alpha = 0.6
  )

  # Clonal vs Subclonal points
  if (subclone$frac1_A == 1) {
    q <- q + ggplot2::geom_point(
      data = subclone,
      ggplot2::aes(
        rlang::.data$nMaj1_A,
        rlang::.data$nMin1_A
      ), size = 5
    )
  } else {
    target_cols <- grep("nM.{5}$|^frac.{3}$", colnames(subclone))
    sol_matrix <- matrix(unlist(subclone[, target_cols]), byrow = TRUE, ncol = 3)
    solutions <- cbind(sol_matrix, rep(1:6, each = 2))
    solutions <- solutions[12:1, ]
    colnames(solutions) <- c("nMaj", "nMin", "frac", "sol")

    solutions_df <- stats::na.omit(as.data.frame(solutions))

    q <- q + ggplot2::geom_point(
      data = solutions_df,
      ggplot2::aes(
        x = rlang::.data$nMaj,
        y = rlang::.data$nMin,
        size = rlang::.data$frac,
        colour = factor(rlang::.data$sol)
      ),
      alpha = 0.75,
      position = ggplot2::position_jitter(width = .05, height = .05),
      shape = 79
    ) +
      ggplot2::scale_size_continuous(guide = "none", limits = c(0, 1), range = c(2, 10)) +
      ggplot2::scale_color_discrete(name = "solution")
  }

  # Final markers
  q <- q + ggplot2::geom_point(ggplot2::aes(x = rlang::.data$nMajcalc, y = rlang::.data$nMincalc), size = 4, shape = 88)
  q <- q + ggplot2::labs(title = paste0(tumourname, " chr", subclone$chr, ": ", subclone$startpos, "-", subclone$endpos))

  print(q)
  grDevices::dev.off()
}

#' Smooth data by running median
#'
#' @param chromosome Denominator on which chromosome each data point belongs. Smoothing is done separately per chromosome
#' @param data The to be smoothed data vector
#' @param k The size of window to be used to take the median over
#' @return A single vector with the smoothed data
#' @author sd11
#' @noRd
runmed_data <- function(chromosome, data, k = 101) {
  data_smoothed <- rep(NA, length(data))
  for (chrom in unique(chromosome)) {
    data_smoothed[chromosome == chrom] <- stats::runmed(data[chromosome == chrom], k)
  }
  return(data_smoothed)
}

#' Plot total copy number split per chromosome
#'
#' This plot contains estimated total copy number from logR, the copy number fit in different colours and a few general stats.
#' It is meant as a single figure replacement for the per chromosome subclones.png figures that can be used for refitting.
#' @param samplename Name of the sample for the plot title
#' @param subclones A subclones.txt file read in as a data.frame
#' @param logr The raw logR read in as a data.frame
#' @param outputfile Full path of file where the figure is to be stored
#' @param purity The samples purity estimate
#' @author sd11
#' @export
totalcn_chrom_plot <- function(
  samplename,
  subclones,
  logr,
  outputfile,
  purity
) {
  # Using data.table::setnames to avoid copying the whole table
  data.table::setnames(logr, 3, "raw_logr")

  # collapse::fcompute/fmutate is faster for smoothing across chromosomes
  # Assuming runmed_data is your custom function
  logr$logr_smoothed <- runmed_data(logr$Chromosome, logr$raw_logr, 101)

  subclones$len <- (subclones$endpos - subclones$startpos) / 1000
  subclones$total_major <- calc_total_cn_major(subclones)
  subclones$total_minor <- calc_total_cn_minor(subclones)
  subclones$total_cn <- subclones$total_minor + subclones$total_major
  subclones$is_subclonal <- subclones$frac1_A < 1
  subclones$is_50_50 <- subclones$frac1_A >= 0.48 & subclones$frac1_A <= 0.52

  ploidy <- calc_ploidy(subclones)
  psi <- psit2psi(purity, ploidy)

  # Convert to data.table if they aren't already
  data.table::setDT(logr)
  data.table::setDT(subclones)
  logr$Position_end <- logr$Position

  # Calculate segment constants once (Vectorized)
  subclones$target_total_cn <- purity * calculate_bb_total_cn(subclones) + 2 * (1 - purity)

  # Set keys for foverlaps (Standard requirement for range joins)
  data.table::setkeyv(subclones, c("chr", "startpos", "endpos"))

  # Perform the join - This maps the correct 'target_total_cn' to every SNP
  logr_joined <- data.table::foverlaps(
    logr,
    subclones,
    by.x = c("Chromosome", "Position", "Position_end"),
    by.y = c("chr", "startpos", "endpos"),
    type = "within",
    nomatch = NA
  )

  # Vectorized calculation of CN columns on the joined data
  logr_joined$total_cn <- logr2tumcn(purity, logr_joined$target_total_cn, logr_joined$logr_smoothed)
  logr_joined$total_cn_psi <- logr2tumcn(purity, psi, logr_joined$logr_smoothed)

  # Replace .N with standard nrow() indexing
  sample_idx <- seq(from = 1, to = nrow(logr_joined), by = 100)
  logr_plot <- logr_joined[sample_idx, ]

  # mixedsort handles the chr1, chr2, chr10 order correctly
  chr_levels <- gtools::mixedsort(unique(as.character(logr_plot$Chromosome)))
  logr_plot$Chromosome <- factor(logr_plot$Chromosome, levels = chr_levels)
  subclones$Chromosome <- factor(subclones$chr, levels = chr_levels)

  max_cn_plot_data <- ceiling(
    collapse::fquantile(
      logr_plot$total_cn_psi, 0.98,
      na.rm = TRUE
    )
  )

  # Optimization: Use weighted quantile for the fit instead of rep() + unlist()
  # This saves massive amounts of memory
  max_cn_plot_fit <- ceiling(
    collapse::fquantile(
      subclones$total_cn, 0.98,
      w = subclones$len, na.rm = TRUE
    )
  )

  max_cn_plot <- max(4, max_cn_plot_data, max_cn_plot_fit, na.rm = TRUE)
  maxpos <- max(logr$Position)

  # 7. Background Data Preparation
  bg_y <- seq(0, max_cn_plot, 2)
  background <- data.frame(
    xmin = 0,
    xmax = maxpos,
    ymin = bg_y + 0.5,
    ymax = bg_y + 1.5
  )

  # 8. Plot Annotations
  prop_subclonal <- round(
    sum(subclones$len[subclones$is_subclonal]) / sum(subclones$len), 2
  )
  homdel <- sum(subclones$len[subclones$total_cn == 0] / 1000)

  plot_subtitle <- paste0(
    "Purity: ", round(purity, 2),
    " - Ploidy: ", round(ploidy, 2),
    " - Hom del: ", round(homdel, 2), "Mb",
    " - Prop. subclonal: ", prop_subclonal
  )
  rect_height_padding <- 0.2

  # Build the actual plot - CNA segments are drawn separately depending on their category as categories have different colours
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = background,
      ggplot2::aes(
        xmin = rlang::.data$xmin,
        xmax = rlang::.data$xmax,
        ymin = rlang::.data$ymin,
        ymax = rlang::.data$ymax
      ),
      fill = "gray80", alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = logr_plot,
      mapping = ggplot2::aes(
        x = rlang::.data$Position,
        y = rlang::.data$total_cn_psi
      ),
      size = 0.5
    ) +
    ggplot2::ylab("Copy Number") +
    ggplot2::scale_y_continuous(breaks = seq(0, max_cn_plot, 2)) +
    # Axis ticks every 10Mb
    ggplot2::scale_x_continuous(
      breaks = seq(1, max(logr$Position), 10000000)[-1],
      labels = round(seq(0, maxpos, 10000000) / 1000000)[-1], expand = c(0, 0)
    ) +
    # Don't restrict the plotting area, zoom. that way segments that go outside the limits are partially plotted still
    ggplot2::coord_cartesian(
      ylim = c(-rect_height_padding, max_cn_plot + rect_height_padding)
    ) +
    ggplot2::facet_wrap(~ rlang::.data$Chromosome, ncol = 2, strip.position = "right") +
    ggplot2::ggtitle(
      bquote(
        atop(
          .(samplename),
          atop(.(plot_subtitle), "")
        )
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        colour = "black", size = 16, face = "plain"
      ),
      axis.text.y = ggplot2::element_text(
        colour = "black", size = 16, face = "plain"
      ),
      axis.title.y = ggplot2::element_text(
        colour = "black", size = 20, face = "plain"
      ),
      strip.text.y = ggplot2::element_text(
        colour = "black", size = 20, face = "plain"
      ),
      plot.title = ggplot2::element_text(
        colour = "black", size = 36, face = "plain", hjust = 0.5
      )
    )

  # Plot the copy number segments - some of the data.frames may be empty, so check for that first before adding to the plot
  sel <- !subclones$is_subclonal
  if (any(sel)) {
    # Minor allele - Normal clonal copy number
    p <- p + ggplot2::geom_rect(
      data = subclones[sel, ],
      mapping = ggplot2::aes(
        xmin = rlang::.data$startpos,
        xmax = rlang::.data$endpos,
        ymin = rlang::.data$total_minor - rect_height_padding,
        ymax = rlang::.data$total_minor + rect_height_padding
      ), fill = "#2f4f4f"
    )
  }
  sel <- subclones$is_subclonal & !subclones$is_50_50
  if (any(sel)) {
    # Minor allele - Normal subclonal copy number
    p <- p + ggplot2::geom_rect(
      data = subclones[sel, ],
      mapping = ggplot2::aes(
        xmin = rlang::.data$startpos,
        xmax = rlang::.data$endpos,
        ymin = rlang::.data$total_minor - rect_height_padding,
        ymax = rlang::.data$total_minor + rect_height_padding
      ), fill = "#2f3f4f"
    )
  }
  sel <- subclones$is_subclonal & subclones$is_50_50
  if (any(sel)) {
    # Minor allele - Subclonal segments right in between two clonal states
    p <- p + ggplot2::geom_rect(
      data = subclones[sel, ],
      mapping = ggplot2::aes(
        xmin = rlang::.data$startpos,
        xmax = rlang::.data$endpos,
        ymin = rlang::.data$total_minor - rect_height_padding,
        ymax = rlang::.data$total_minor + rect_height_padding
      ), fill = "#2f3f4f", colour = "red"
    )
  }
  sel <- !subclones$is_subclonal
  if (any(sel)) {
    # Major allele - clonal copy number
    p <- p + ggplot2::geom_rect(
      data = subclones[sel, ],
      mapping = ggplot2::aes(
        xmin = rlang::.data$startpos,
        xmax = rlang::.data$endpos,
        ymin = rlang::.data$total_cn - rect_height_padding,
        ymax = rlang::.data$total_cn + rect_height_padding
      ), fill = "#E69F00"
    )
  }
  sel <- subclones$is_subclonal
  if (any(sel)) {
    # Major allele - subclonal copy number
    p <- p + ggplot2::geom_rect(
      data = subclones[sel, ],
      mapping = ggplot2::aes(
        xmin = rlang::.data$startpos,
        xmax = rlang::.data$endpos,
        ymin = rlang::.data$total_cn - rect_height_padding,
        ymax = rlang::.data$total_cn + rect_height_padding
      ), fill = "#E55300"
    )
  }

  grDevices::png(outputfile, width = 2000, height = 1300, type = "cairo")
  print(p)
  grDevices::dev.off()
}

#' Plot allele ratios from raw segmented data
#'
#' @param samplename Name of the sample for the plot title
#' @param bafsegmented The BAFsegmented data read in as a data.frame
#' @param logrsegmented The logRsegmented data read in as a data.frame
#' @param outputfile Full path of file where the figure is to be stored
#' @param logr A data.frame with the logr, either logr or allelecounts must be supplied
#' @param max.plot.cn Maximum y-axis value to plot (Default: 5)
#' @author sd11
#' @export
allele_ratio_plot <- function(
  samplename, bafsegmented,
  logrsegmented, outputfile,
  logr, max.plot.cn = 5
) {
  if (nrow(logr) < 2000000) {
    platform <- "SNP6"
  } else {
    platform <- "WGS"
  }

  bafsegmented$Chromosome <- factor(bafsegmented$Chromosome, levels = gtools::mixedsort(unique(bafsegmented$Chromosome)))
  colnames(logrsegmented) <- c("Chromosome", "Position", "logRseg")
  logrsegmented$Chromosome <- factor(logrsegmented$Chromosome, levels = S4Vectors::levels(bafsegmented$Chromosome))

  colnames(logr)[3] <- "raw_logr"
  logr$copy_ratio_binned <- runmed_data(logr$Chromosome, exp(logr$raw_logr))
  logr$Chromosome <- factor(logr$Chromosome, levels = S4Vectors::levels(bafsegmented$Chromosome))
  allelecounts <- logr

  copyratio_binnedLogR <- as.data.frame(array(NA, c(nrow(bafsegmented), 8)))
  colnames(copyratio_binnedLogR) <- c("Chromosome", "Position", "ratioBAF", "ratioBAFphased", "ratioBAF_alt", "ratioBAFphased_alt", "ratioBAFseg", "ratioBAFseg_alt")
  copyratio_binnedLogR$Chromosome <- bafsegmented$Chromosome
  copyratio_binnedLogR$Position <- bafsegmented$Position

  log_info("Calculating copy ratios..")
  for (chrom in unique(bafsegmented$Chromosome)) {
    print(chrom)

    baf_chrom <- bafsegmented[bafsegmented$Chromosome == chrom, ]
    logrseg_chrom <- logrsegmented[logrsegmented$Chromosome == chrom, ]

    baf_sel <- baf_chrom$Position %in% intersect(baf_chrom$Position, logrseg_chrom$Position)
    logrseg_sel <- logrseg_chrom$Position %in% intersect(baf_chrom$Position, logrseg_chrom$Position)
    ratio_sel <- which(copyratio_binnedLogR$Chromosome == chrom)[baf_sel]

    copyratio_binnedLogR$ratioBAFseg[ratio_sel] <- (baf_chrom$BAFseg[baf_sel] * (2^logrseg_chrom$logRseg[logrseg_sel]))
    copyratio_binnedLogR$ratioBAFseg_alt[ratio_sel] <- (-(baf_chrom$BAFseg[baf_sel] - 1) * (2^logrseg_chrom$logRseg[logrseg_sel]))
  }

  background <- data.frame(y = seq(0, max.plot.cn, 0.5))

  log_info("Plotting..")
  if (platform == "WGS") {
    sel <- seq(1, nrow(allelecounts), 100)
  } else {
    sel <- rep(TRUE, nrow(allelecounts))
  }

  plot_title <- samplename
  copy_ratio <- ggplot2::ggplot(allelecounts[sel, ]) +
    ggplot2::geom_hline(
      data = background, mapping = ggplot2::aes(yintercept = rlang::.data$y),
      colour = "black", alpha = 0.3
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = rlang::.data$Position,
        y = rlang::.data$copy_ratio_binned
      ),
      alpha = 0.5, size = 0.9, colour = "darkgreen"
    ) +
    ggplot2::facet_grid(~ rlang::.data$Chromosome, scales = "free_x", space = "free_x")
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::ylim(0, max.plot.cn) +
    ggplot2::ylab("Copy Ratio") +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        colour = "black", size = 18, face = "plain"
      ),
      axis.title.y = ggplot2::element_text(
        colour = "black", size = 20, face = "plain"
      ),
      strip.text.x = ggplot2::element_text(
        colour = "black", size = 16, face = "plain"
      ),
      plot.title = ggplot2::element_text(
        colour = "black", size = 36, face = "plain", hjust = 0.5
      )
    )

  if (platform == "WGS") {
    sel <- seq(1, nrow(copyratio_binnedLogR), 100)
  } else {
    sel <- rep(TRUE, nrow(copyratio_binnedLogR))
  }
  as_copy_ratio_seg <- ggplot2::ggplot(copyratio_binnedLogR[sel, ]) +
    ggplot2::geom_hline(
      data = background,
      mapping = ggplot2::aes(yintercept = rlang::.data$y),
      colour = "black", alpha = 0.3
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = rlang::.data$Position, y = rlang::.data$ratioBAFseg_alt
      ), alpha = 0.5, size = 0.9, colour = "darkblue"
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = rlang::.data$Position, y = rlang::.data$ratioBAFseg
      ), alpha = 0.5, size = 0.9, colour = "purple"
    ) +
    ggplot2::facet_grid(~ rlang::.data$Chromosome, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::ylim(0, max.plot.cn) +
    ggplot2::ylab("AS Copy Ratio - Segm") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        colour = "black", size = 18, face = "plain"
      ),
      axis.title.y = ggplot2::element_text(
        colour = "black", size = 20, face = "plain"
      ),
      strip.text.x = ggplot2::element_text(
        colour = "black", size = 16, face = "plain"
      ),
      plot.title = ggplot2::element_text(
        colour = "black", size = 36, face = "plain"
      )
    )
  grDevices::png(outputfile, width = 2000, height = 750, type = "cairo")
  gridExtra::grid.arrange(
    gridExtra::arrangeGrob(copy_ratio, as_copy_ratio_seg, ncol = 1)
  )
  grDevices::dev.off()
}

#' Plot relative coverage of tumour and normal
#'
#' @param samplename Name of the sample for the plot title
#' @param allelecounts Combined allele counts of tumour and normal, read in as a data.frame
#' @param outputfile Full path of file where the figure is to be stored
#' @param max.y The max Y-axis value to be plotted
#' @author sd11
#' @export
coverage_plot <- function(samplename, allelecounts, outputfile, max.y = 4) {
  log_info("Normalising allele counts..")
  allelecounts$tumour <- allelecounts$mutCountT1 + allelecounts$mutCountT2
  allelecounts$tumour <- allelecounts$tumour / collapse::fmedian(allelecounts$tumour, na.rm = TRUE)
  allelecounts$normal <- allelecounts$mutCountN1 + allelecounts$mutCountN2
  allelecounts$normal <- allelecounts$normal / collapse::fmedian(allelecounts$normal, na.rm = TRUE)

  log_info("Smoothing data..")
  # res = bin_coverage_tumour(allelecounts, binsize=10000)
  # allelecounts$tumour_binned = res$tumour_binned
  allelecounts$tumour_binned <- runmed_data(allelecounts$Chromosome, allelecounts$tumour)

  # res = bin_coverage_normal(allelecounts, binsize=10000)
  # allelecounts$normal_binned = res$normal_binned
  # rm(res)
  allelecounts$normal_binned <- runmed_data(allelecounts$Chromosome, allelecounts$normal)
  allelecounts$Chromosome <- factor(allelecounts$Chromosome, levels = gtools::mixedsort(unique(allelecounts$Chromosome)))

  background <- data.frame(y = seq(0, 2, 0.5))
  plot_title <- samplename
  p <- ggplot2::ggplot(allelecounts[seq(1, nrow(allelecounts), 100), ]) +
    ggplot2::geom_hline(
      data = background, mapping = ggplot2::aes(yintercept = rlang::.data$y),
      colour = "black", alpha = 0.3
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = rlang::.data$Position, y = rlang::.data$normal_binned),
      alpha = 0.5, size = 0.5, colour = "darkgreen"
    ) +
    ggplot2::facet_grid(~Chromosome, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::ylab("Normal") +
    ggplot2::scale_y_continuous(breaks = c(0:2), limits = c(0, 2)) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        colour = "black", size = 18, face = "plain"
      ),
      axis.title.y = ggplot2::element_text(
        colour = "black", size = 20, face = "plain"
      ),
      strip.text.x = ggplot2::element_text(
        colour = "black", size = 16, face = "plain"
      ),
      plot.title = ggplot2::element_text(
        colour = "black", size = 36, face = "plain", hjust = 0.5
      )
    )

  background <- data.frame(y = seq(0, max.y, 0.5))
  p3 <- ggplot2::ggplot(allelecounts[seq(1, nrow(allelecounts), 100), ]) +
    ggplot2::geom_hline(
      data = background,
      mapping = ggplot2::aes(yintercept = rlang::.data$y),
      colour = "black", alpha = 0.3
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = rlang::.data$Position,
        y = rlang::.data$tumour_binned
      ),
      alpha = 0.5, size = 0.5, colour = "darkgreen"
    ) +
    ggplot2::facet_grid(~Chromosome, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::ylim(0, max.y) +
    ggplot2::ylab("Tumour") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        colour = "black", size = 18, face = "plain"
      ),
      axis.title.y = ggplot2::element_text(
        colour = "black", size = 20, face = "plain"
      ),
      strip.text.x = ggplot2::element_text(
        colour = "black", size = 16, face = "plain"
      ),
      plot.title = ggplot2::element_text(
        colour = "black", size = 36, face = "plain"
      )
    )
  grDevices::png(outputfile, width = 2000, height = 750, type = "cairo")
  gridExtra::grid.arrange(gridExtra::arrangeGrob(p, p3, ncol = 1))
  grDevices::dev.off()
}
