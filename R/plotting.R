#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
create.haplotype.plot = function(chrom.position, points.blue, points.red, x.min, x.max, title, xlab, ylab) {
  par(pch=".", cex=1, cex.main=0.8, cex.axis = 0.6, cex.lab=0.7,yaxp=c(-0.05,1.05,6))
  plot(c(x.min,x.max), c(0,1), type="n", main=title, xlab=xlab, ylab=ylab)
  if (length(chrom.position) > 0) {
    points(chrom.position, points.blue, col="blue")
    points(chrom.position, points.red, col="red")
  }
}

#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
create.segmented.plot = function(chrom.position, points.red, points.green, x.min, x.max, title, xlab, ylab, prior_bkps_pos=NULL) {
  par(mar = c(5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, cex.lab = 2)
  plot(c(x.min,x.max), c(0,1), pch=".", type="n", main=title, xlab=xlab, ylab=ylab)
  points(chrom.position, points.red, pch=".", col="red", cex=2)
  points(chrom.position, points.green, pch=19, cex=0.5, col="green")
  if (!is.null(prior_bkps_pos)) {
    for (i in 1:length(prior_bkps_pos)) {
      abline(v=prior_bkps_pos[i])
    }
  }
}

#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
create.baf.plot = function(chrom.position, points.red.blue, plot.red, points.darkred, points.darkblue, x.min, x.max, title, xlab, ylab, prior_bkps_pos=NULL) {
  par(mar = c(5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, cex.lab = 2)
  plot(c(x.min,x.max), c(0,1), pch=".", type = "n", main=title, xlab=xlab, ylab=ylab)
  points(chrom.position, points.red.blue, pch=".", col=ifelse(plot.red, "red", "blue"), cex=2)
  points(chrom.position, points.darkred, pch=19, cex=0.5, col="darkred")
  points(chrom.position, points.darkblue, pch=19, cex=0.5, col="darkblue")
  if (!is.null(prior_bkps_pos)) {
    for (i in 1:length(prior_bkps_pos)) {
      abline(v=prior_bkps_pos[i])
    }
  }
}

#' Function that creates the plots for subclonal copy number
#' Note: This is a plot PER chromosome.
#' @noRd
create.subclonal.cn.plot = function(chrom, chrom.position, LogRposke, LogRchr, BAFchr, BAFsegchr, BAFpvalschr, subcloneres, siglevel, x.min, x.max, title, xlab, ylab.logr, ylab.baf, breakpoints_pos=NULL, svs_pos=NULL) {

  plot_breakpoints = function(breakpoints, svs_pos) {
    # Plot the breakpoints
    if (!is.null(breakpoints)) {
      for (i in 1:length(breakpoints)) {
        abline(v=breakpoints[i], col="darkgrey", lwd=1)
      }
    }
    
    # Overplot the SV breakpoints, if supplied
    if (!is.null(svs_pos)) {
      for (i in 1:length(svs_pos)) {
        abline(v=svs_pos[i], lty=3, col="lightgreen", lwd=1)
      }
    }
  }
  
  # Plot the logR
  par(mar=c(2.5,2.5,2.5,0.25), cex=0.4, cex.main=1.5, cex.axis=1, cex.lab=1, mfrow=c(2,1))
  plot(c(x.min, x.max), c(-3,3), pch=".", type="n", main=title, xlab=xlab, ylab=ylab.logr)
  points(LogRposke/1000000, LogRchr, pch=".", col="grey")
  plot_breakpoints(breakpoints_pos, svs_pos)
  
  # Plot BAF
  plot(c(x.min, x.max), c(0,1), pch=".", type="n", main=title, xlab=xlab, ylab=ylab.baf)
  points(chrom.position, BAFchr, pch=".", col="grey")
  plot_breakpoints(breakpoints_pos, svs_pos)
  
  # Plot segments in top of BAF
  points(chrom.position, BAFsegchr, pch=19, cex=0.5, col=ifelse(BAFpvalschr>siglevel, "darkgreen", "red"))
  points(chrom.position, 1-BAFsegchr, pch=19, cex=0.5, col=ifelse(BAFpvalschr>siglevel, "darkgreen", "red"))
  for (i in 1:dim(subcloneres)[1]) {
    if(subcloneres[i,1]==chrom) {
      text((as.numeric(subcloneres[i,"startpos"])+as.numeric(subcloneres[i,"endpos"]))/2/1000000,as.numeric(subcloneres[i,"BAF"])-0.04,
           paste(subcloneres[i,"nMaj1_A"],"+",subcloneres[i,"nMin1_A"],": ",100*round(as.numeric(subcloneres[i,"frac1_A"]),3),"%",sep=""),cex = 0.8)
      if(!is.na(subcloneres[i,"nMaj2_A"])) {
        text((as.numeric(subcloneres[i,"startpos"])+as.numeric(subcloneres[i,"endpos"]))/2/1000000,as.numeric(subcloneres[i,"BAF"])-0.08,
             paste(subcloneres[i,"nMaj2_A"],"+",subcloneres[i,"nMin2_A"],": ",100*round(as.numeric(subcloneres[i,"frac2_A"]),3),"%",sep=""), cex = 0.8)
      }
    }
  }
}


#' Make GW CN plot where subclonal CN is represented as a mixture of two states
#' NAP - July 2020 - updated main title now replacing 'cellularity' with 'purity' and 'goodness-of-fit' with 'PGAclonal' + adding TUMOURNAME
#' @noRd
create.bb.plot.average = function(bafsegmented, ploidy, rho, goodnessOfFit, pos_min, pos_max, segment_states_min, segment_states_tot, chr.segs, chr.names, tumourname, ylim=5) {
  # Plot main frame and title
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  maintitle = paste0(substring(tumourname, 40, first = T),", Ploidy: ",sprintf("%1.2f",ploidy),", Purity: ",sprintf("%2.0f",rho*100),"%, PGAclonal: ",sprintf("%2.1f",goodnessOfFit*100),"%")
  #maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit*100),"%",sep="")
  plot(c(1,nrow(bafsegmented)), c(0,ylim), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
  abline(v=0,lty=1,col="lightgrey")
  # Horizontal lines for y=0 to y=5
  abline(h=c(0:ylim),lty=1,col="lightgrey")
  # Minor allele in gray, total CN in orange
  segments(x0=pos_min, y0=segment_states_min, x1=pos_max, y1=segment_states_min, col="#2f4f4f", pch="|", lwd=6, lend=1)
  segments(x0=pos_min, y0=segment_states_tot, x1=pos_max, y1=segment_states_tot, col="#E69F00", pch="|", lwd=6, lend=1)

  # Plot the vertical lines that show start/end of a chromosome
  chrk_tot_len = 0
  for (i in 1:length(chr.segs)) {
    chrk = chr.segs[[i]];
    chrk_hetero = names(bafsegmented)[chrk]
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,ylim,chr.names[i], pos = 1, cex = 2)
    abline(v=vpos,lty=1,col="lightgrey")
  }
}

#' Make GW CN plot where subclonal CN is represented by two separate states
#' NAP - July 2020 - updated main title now replacing 'cellularity' with 'purity' and 'goodness-of-fit' with 'PGAclonal' + adding TUMOURNAME
#' @noRd
create.bb.plot.subclones = function(bafsegmented, subclones, ploidy, rho, goodnessOfFit, pos_min, pos_max, subcl_min, subcl_max, is_subclonal, is_subclonal_maj, is_subclonal_min, chr.segs, chr.names, tumourname, ylim=5) {
	par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
	maintitle = paste0(substring(tumourname, 40, first = T),", Ploidy: ",sprintf("%1.2f",ploidy),", Purity: ",sprintf("%2.0f",rho*100),"%, PGAclonal: ",sprintf("%2.1f",goodnessOfFit*100),"%")
        # maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit*100),"%",sep="")
	plot(c(1,nrow(bafsegmented)), c(0,ylim), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
	abline(v=0,lty=1,col="lightgrey")
	# Minor allele clonal and lowest of the two states when subclonal
	segments(x0=pos_min, y0=subclones$nMin1_A-0.1, 
		 x1=pos_max, y1=subclones$nMin1_A-0.1, col="#2f4f4f", pch="|", 
		 lwd=ifelse(is_subclonal_min, 6*subclones$frac1_A, 6), lend=1)

	if (sum(is_subclonal) > 0) {
		# Minor allele highest of the two states when subclonal
		segments(x0=subcl_min, y0=subclones$nMin2_A[is_subclonal]-0.1, 
			 x1=subcl_max, y1=subclones$nMin2_A[is_subclonal]-0.1, col="#2f4f4f", pch="|", 
			 lwd=ifelse(is_subclonal_min[is_subclonal], 6*subclones$frac2_A[is_subclonal], 0), lend=1)

		# Total CN, when minor allele subclonal CN (one of the two alleles)
		segments(x0=subcl_min, y0=subclones$nMaj1_A[is_subclonal]+subclones$nMin1_A[is_subclonal]+0.1, 
			 x1=subcl_max, y1=subclones$nMaj1_A[is_subclonal]+subclones$nMin1_A[is_subclonal]+0.1, col="#E69F00", pch="|", 
			 lwd=ifelse(is_subclonal_min[is_subclonal], 6*subclones$frac1_A[is_subclonal], 0), lend=1)

		# Total CN, when minor allele subclonal CN (the other allele)
		segments(x0=subcl_min, y0=subclones$nMaj2_A[is_subclonal]+subclones$nMin2_A[is_subclonal]+0.1, 
			 x1=subcl_max, y1=subclones$nMaj2_A[is_subclonal]+subclones$nMin2_A[is_subclonal]+0.1, col="#E69F00", pch="|", 
			 lwd=ifelse(is_subclonal_min[is_subclonal], 6*subclones$frac2_A[is_subclonal], 0), lend=1)
	}

	# Total CN, when major allele clonal and subclonal, unless the minor allele is subclonal (then plot nothing, done above)
	segments(x0=pos_min, y0=subclones$nMaj1_A+subclones$nMin1_A+0.1, 
		 x1=pos_max, y1=subclones$nMaj1_A+subclones$nMin1_A+0.1, col="#E69F00", pch="|", 
		 lwd=ifelse(is_subclonal_maj & (!is_subclonal_min), 6*subclones$frac1_A, 0), lend=1)

	# Total CN, when subclonal major allele and not subclonal minor allele (the other allele)
	segments(x0=pos_min, y0=subclones$nMaj2_A+subclones$nMin2_A+0.1, 
		 x1=pos_max, y1=subclones$nMaj2_A+subclones$nMin2_A+0.1, col="#E69F00", pch="|", 
		 lwd=ifelse(is_subclonal_maj & (!is_subclonal_min), 6*subclones$frac2_A, 0), lend=1)

	# Total allele when major and minor both non-subclonal
	segments(x0=pos_min, y0=subclones$nMaj1_A+subclones$nMin1_A+0.1, 
		 x1=pos_max, y1=subclones$nMaj1_A+subclones$nMin1_A+0.1, col="#E69F00", pch="|", 
		 lwd=ifelse((!is_subclonal_maj) & (!is_subclonal_min), 6, 0), lend=1)

	chrk_tot_len = 0
	for (i in 1:length(chr.segs)) {
		chrk = chr.segs[[i]];
		chrk_hetero = names(bafsegmented)[chrk]
		chrk_tot_len_prev = chrk_tot_len
		chrk_tot_len = chrk_tot_len + length(chrk_hetero)
		vpos = chrk_tot_len;
		tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
		text(tpos,ylim,chr.names[i], pos = 1, cex = 2)
		abline(v=vpos,lty=1,col="lightgrey")
	}
}

#' Code extracted from the plot in clonal_ascat find_centroid_of_global_minima.
#' Note: This is a temporary function and VERY similar to clonal_runascat.plot1()
#' @noRd
#'
clonal_findcentroid.plot = function(minimise, dist_choice, d, psis, rhos, new_bounds) {
  par(mar = c(5,5,0.5,0.5), cex=0.75, cex.lab=2, cex.axis=2)
  if(minimise){ #DCW 240314 reverse colour palette, so blue always corresponds to best region
    hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
  } else {
    hmcol = colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256)
  }
  if ( dist_choice == 4 ) {
    image(d, col = hmcol, axes = F, xlab = "Ploidy", ylab = "Aberrant cell fraction")
  } else  {
    image(log(d), col = hmcol, axes = F, xlab = "Ploidy", ylab = "Aberrant cell fraction")
  }
  psi_min = new_bounds$psi_min
  psi_max = new_bounds$psi_max
  rho_min = new_bounds$rho_min
  rho_max = new_bounds$rho_max
  
  psi_range = psi_max - psi_min
  rho_range = rho_max - rho_min

  psi_min_label = ceiling( 10 * psi_min )/10
  psi_max_label = floor( 10 * psi_max )/10
  psi_label_interval = 0.1
  
  psi_min_label_standardised = ( psi_min_label - psi_min ) / psi_range
  psi_max_label_standardised = ( psi_max_label - psi_min ) / psi_range
  psi_label_interval_standardised = psi_label_interval / psi_range
  
  rho_min_label = ceiling( 100 * rho_min )/100
  rho_max_label = floor( 100 * rho_max )/100
  rho_label_interval = 0.01
  
  rho_min_label_standardised = ( rho_min_label - rho_min ) / rho_range
  rho_max_label_standardised = ( rho_max_label - rho_min ) / rho_range
  rho_label_interval_standardised = rho_label_interval / rho_range
  
  axis(1, at = seq(psi_min_label_standardised, psi_max_label_standardised, by = psi_label_interval_standardised), labels = seq(psi_min_label, psi_max_label, by = psi_label_interval))
  axis(2, at = seq(rho_min_label_standardised, rho_max_label_standardised, by = rho_label_interval_standardised), labels = seq(rho_min_label, rho_max_label, by = rho_label_interval))
  
  points( ( psis - psi_min ) / psi_range , ( rhos - rho_min ) / rho_range , col=c("green", "darkgreen"), pch="X", cex = 2 )
}


#' Plot Battenberg copy number solutions for a segment 
#' 
#' \code{squaresplot} plots the different Battenberg copy number solutions for a segment
#' 
#' The plot is output to the run directory as "tumourname_squares_chr_position.png/pdf"
#' 
#' @param tumourname Sample name
#' @param run_dir Running directory
#' @param segment_chr Chromosome containing the segment to be investigated
#' @param segment_pos Chromosomal position within the segment in Mb (e.g. 90M)
#' @param platform_gamma Platform-specific gamma value (0.55 for SNP6, 1 for NGS), default 1
#' @param pdf Output format: 0 for png (default), 1 for pdf
#' @param binwidth_baf BAF isobafline spacing, default 0.25
#' @param xylimits x/y-axis limits, default c(-0.2,5)
#' @author jd
#' @export
squaresplot <- function(tumourname, run_dir, segment_chr, segment_pos, platform_gamma=1, pdf=0, binwidth_baf=0.25, xylimits=c(-0.2,5)) {
  
  if (pdf)
    pdf(file = paste(run_dir,tumourname,"_squares","_chr",segment_chr,"_",segment_pos,".pdf", sep=""), width = 7, height = 7)
  else
    png(filename = paste(run_dir,tumourname,"_squares","_chr",segment_chr,"_",segment_pos,".png", sep=""), width = 1200, height = 1200, res = 200)
  
  # read in and augment data
  segment_pos <- as.numeric(gsub("M", "000000", segment_pos))
  subclones <- read.table(paste(run_dir, tumourname, "_subclones.txt", sep=""), header=T, stringsAsFactors=F)
  subclone <- subclones[(subclones$chr == segment_chr) & (subclones$startpos <= segment_pos) & (subclones$endpos >= segment_pos),]
  rhopsi <- read.table(paste(run_dir, tumourname, "_rho_and_psi.txt", sep=""), header = T, stringsAsFactors=F)
  rhopsi <- rhopsi[which(rhopsi$is.best == TRUE), c("rho", "psi")] 
  
  nMincalc <- (rhopsi$rho-1-(subclone$BAF-1)*2^(subclone$LogR/platform_gamma)*((1-rhopsi$rho)*2+rhopsi$rho*rhopsi$psi))/rhopsi$rho
  nMajcalc <- (rhopsi$rho-1+subclone$BAF*2^(subclone$LogR/platform_gamma)*((1-rhopsi$rho)*2+rhopsi$rho*rhopsi$psi))/rhopsi$rho

  subclone <- data.frame(subclone, rhopsi, nMincalc, nMajcalc)
  
  # helper function to calculate isobaflines
  isobafline <- function(nB, cstbaf) {
    (1-rhopsi$rho+rhopsi$rho*nB-cstbaf*(2-2*rhopsi$rho)-rhopsi$rho*cstbaf*nB)/(rhopsi$rho*cstbaf)
  }
  
  # create grid for allelic copynumber
  ngrid <- data.frame(nMaj=seq(0,5,1), nMin=seq(0,5,1))
  
  # start plotting - setup
  q <- ggplot2::ggplot(data = ngrid, aes(nMaj, nMin)) + ggplot2::scale_x_continuous(breaks=0:max(xylimits), limits=xylimits) + ggplot2::scale_y_continuous(breaks=0:max(xylimits), limits=xylimits) + ggplot2::coord_fixed()
  q <- q + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_line(colour="darkgrey", size = 0.5), panel.grid.minor = ggplot2::element_blank())
  
  # add isobaflines
  for (bafval in seq(0,1,binwidth_baf)) {
    q <- q + ggplot2::stat_function(fun = isobafline, args = list(cstbaf = bafval), colour="blue", alpha=0.6)
  }
  q <- q + ggplot2::stat_function(fun = isobafline, args = list(cstbaf = subclone$BAF), colour="green")
  
  # add isologrline
  df = data.frame(flnMaj = floor(nMajcalc)-0.2, cnMin = ceiling(nMincalc)+0.2, cnMaj = ceiling(nMajcalc)+0.2, flnMin = floor(nMincalc)-0.2)
  q <- q + ggplot2::geom_segment(data = df,
                        aes(x = flnMaj, y = cnMin, xend = cnMaj, yend = flnMin), colour="red", alpha = 0.6)
  
  # if clonal segment, only plot clonal solution
  if (subclone$frac1_A == 1) {
    q <- q + ggplot2::geom_point(data = subclone, aes(nMaj1_A, nMin1_A), size = 5)
  } else { # if subclonal, plot all equivalent solutions
    solutions <- matrix(unlist(subclone[,grep("nM.{5}$|^frac.{3}$", colnames(subclone))]), byrow = T, ncol = 3)
    solutions <- cbind(solutions, rep(1:6,rep(2,6)))[12:1,]
    colnames(solutions) <- c("nMaj", "nMin", "frac", "sol")
    solutions <- na.omit(as.data.frame(solutions))
    q <- q + ggplot2::geom_point(data = solutions, aes(nMaj, nMin, size=frac, colour=factor(sol)), alpha=0.75, position = ggplot2::position_jitter(width = .05, height = .05), shape = 79) +
      ggplot2::scale_size_continuous(guide=F, limits=c(0,1) ,range = c(2,10)) + ggplot2::scale_color_discrete(name="solution")
  }
  
  # plot precise values, as calculated by battenberg
  q <- q + ggplot2::geom_point(data=subclone, aes(nMajcalc, nMincalc), size=4, shape = 88)
  q <- q + ggplot2::labs(title = paste(tumourname," chr",subclone$chr,": ",subclone$startpos,"-",subclone$endpos, sep=""))
  print(q)
  dev.off()
}

#' Smooth data by running median
#' 
#' @param chromosome Denominator on which chromosome each data point belongs. Smoothing is done separately per chromosome
#' @param data The to be smoothed data vector
#' @param k The size of window to be used to take the median over
#' @return A single vector with the smoothed data
#' @author sd11
#' @noRd
runmed_data = function(chromosome, data, k=101) {
  data_smoothed = rep(NA, length(data))
  for (chrom in unique(chromosome)) {
   data_smoothed[chromosome==chrom] = runmed(data[chromosome==chrom], k)
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
totalcn_chrom_plot = function(samplename, subclones, logr, outputfile, purity) {
  
  # Smooth the logR
  colnames(logr)[3] = "raw_logr"
  logr$logr_smoothed = runmed_data(logr$Chromosome, logr$raw_logr, 101)
  
  # Prepare subclones data
  subclones$len = subclones$endpos/1000-subclones$startpos/1000
  subclones$total_major = calc_total_cn_major(subclones)
  subclones$total_minor = calc_total_cn_minor(subclones)
  subclones$total_cn = subclones$total_minor + subclones$total_major
  subclones$is_subclonal = subclones$frac1_A < 1
  subclones$is_50_50 = subclones$frac1_A >= 0.48 & subclones$frac1_A <= 0.52
  
  # Calculate psi from the data
  ploidy = calc_ploidy(subclones)
  psi = psit2psi(purity, ploidy)

  # Estimate total CN for each segment based on the logR
  logr$total_cn = NA
  logr$total_cn_psi = NA
  for (i in (1:nrow(subclones))) {
    print(i)
    sel = which(logr$Chromosome == subclones$chr[i] & logr$Position >= subclones$startpos[i] & logr$Position <= subclones$endpos[i])
    tumour_cn = calculate_bb_total_cn(subclones[i,,drop=F])
    total_cn = purity*tumour_cn + 2*(1-purity)
    logr$total_cn[sel] = logr2tumcn(purity, total_cn, logr$logr_smoothed[sel])
    logr$total_cn_psi[sel] = logr2tumcn(purity, psi, logr$logr_smoothed[sel])
  }
  
  # Plot every 100 data point, there are too many for them all to be seen
  logr_plot = logr[seq(1, nrow(logr), 100),]
  
  # Sync the levels for chromosome so that all corresponding data ends up in the same plot
  logr_plot$Chromosome = factor(logr_plot$Chromosome, levels=gtools::mixedsort(unique(logr_plot$Chromosome)))
  subclones$Chromosome = factor(subclones$chr, levels=levels(logr_plot$Chromosome))
  
  # Set plot boundaries for x and y - take as y value the maximum between the data and the fit
  max_cn_plot_data = ceiling(quantile(logr_plot$total_cn_psi, c(.98), na.rm=T))
  max_cn_plot_fit = ceiling(quantile(unlist(lapply(1:nrow(subclones), function(i) rep(subclones$total_cn[i], subclones$len[i]))), c(.98), na.rm=T))
  max_cn_plot = ifelse(max_cn_plot_fit > max_cn_plot_data, max_cn_plot_fit, max_cn_plot_data)
  maxpos = max(logr$Position)
  
  # catch case when there is no clonal CNA called
  if (is.na(max_cn_plot) | max_cn_plot < 4) {
    max_cn_plot = 4
  }
  
  # These are the grey lines in the background
  background = data.frame(xmin=rep(0, (max_cn_plot/2)+1), 
                          xmax=rep(max(logr$Position), (max_cn_plot/2)+1),
                          ymin=seq(0, max_cn_plot, 2)+0.5, 
                          ymax=seq(0, max_cn_plot, 2)+1.5)

  # Calc a couple of stats for the plot title
  genome_50_50 = sum(subclones$len[subclones$is_50_50]/1000)
  prop_subclonal = round(sum(subclones$len[subclones$is_subclonal]) / sum(subclones$len), 2)
  homdel = sum(subclones$len[subclones$total_cn == 0]/1000)
  plot_title = samplename
  plot_subtitle = paste0("Purity: ", round(purity, 2), " - Ploidy: ", round(ploidy, 2), " - Hom del: ", round(homdel, 2), "Mb - Prop. subclonal: ", prop_subclonal, " - Subclonal 50/50: ", round(genome_50_50, 2), "Mb")
  
  rect_height_padding = 0.2
  
  # Build the actual plot - CNA segments are drawn separately depending on their category as categories have different colours
  p = ggplot() + 
    geom_rect(data=background, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.5) +
    geom_point(data=logr_plot, mapping=aes(x=Position, y=total_cn_psi), size=0.5) + 
    ylab("Copy Number") +
    scale_y_continuous(breaks=seq(0, max_cn_plot, 2)) +  #, limits=c(-rect_height_padding, max_cn_plot+rect_height_padding)
    # Axis ticks every 10Mb
    scale_x_continuous(breaks=seq(1, max(logr$Position), 10000000)[-1], labels=round(seq(0, maxpos, 10000000) / 1000000)[-1], expand=c(0, 0)) +
    # Don't restrict the plotting area, zoom. that way segments that go outside the limits are partially plotted still
    coord_cartesian(ylim=c(-rect_height_padding, max_cn_plot+rect_height_padding)) +
    facet_wrap(~Chromosome, ncol=2, strip.position="right") + 
    # ggtitle(plot_title) +
    ggtitle(bquote(atop(.(plot_title), atop(.(plot_subtitle), "")))) +
    theme_bw() + theme(axis.title.x=element_blank(),
                       axis.text.x=element_text(colour="black",size=16,face="plain"),
                       axis.text.y = element_text(colour="black",size=16,face="plain"), 
                       axis.title.y = element_text(colour="black",size=20,face="plain"),
                       strip.text.y = element_text(colour="black",size=20,face="plain"),
                       plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5))
  
  # Plot the copy number segments - some of the data.frames may be empty, so check for that first before adding to the plot
  sel = !subclones$is_subclonal
  if (any(sel)) {
    # Minor allele - Normal clonal copy number
    p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_minor-rect_height_padding, ymax=total_minor+rect_height_padding), fill="#2f4f4f")
  }
  sel = subclones$is_subclonal & !subclones$is_50_50
  if (any(sel)) {
    # Minor allele - Normal subclonal copy number
    p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_minor-rect_height_padding, ymax=total_minor+rect_height_padding), fill="#2f3f4f")
  }
  sel = subclones$is_subclonal & subclones$is_50_50
  if (any(sel)) {
    # Minor allele - Subclonal segments right in between two clonal states
    p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_minor-rect_height_padding, ymax=total_minor+rect_height_padding), fill="#2f3f4f", colour="red")
  }
  sel = !subclones$is_subclonal
  if (any(sel)) {
    # Major allele - clonal copy number
    p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_cn-rect_height_padding, ymax=total_cn+rect_height_padding), fill="#E69F00")
  }
  sel = subclones$is_subclonal
  if (any(sel)) {
    # Major allele - subclonal copy number
    p = p + geom_rect(data=subclones[sel, ], mapping=aes(xmin=startpos, xmax=endpos, ymin=total_cn-rect_height_padding, ymax=total_cn+rect_height_padding), fill="#E55300")
  }
  
  png(outputfile, width=2000, height=1300)
  print(p)
  dev.off()
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
allele_ratio_plot = function(samplename, bafsegmented, logrsegmented, outputfile, logr, max.plot.cn=5) {

  if (nrow(logr) < 2000000) {
	  platform = "SNP6"
  } else {
	  platform = "WGS"
  }

  bafsegmented$Chromosome = factor(bafsegmented$Chromosome, levels=gtools::mixedsort(unique(bafsegmented$Chromosome)))
  colnames(logrsegmented) = c("Chromosome", "Position", "logRseg")
  logrsegmented$Chromosome = factor(logrsegmented$Chromosome, levels=levels(bafsegmented$Chromosome))
  
  colnames(logr)[3] = "raw_logr"
  logr$copy_ratio_binned = runmed_data(logr$Chromosome, exp(logr$raw_logr))
  logr$Chromosome = factor(logr$Chromosome, levels=levels(bafsegmented$Chromosome))
  allelecounts = logr

  copyratio_binnedLogR = as.data.frame(array(NA, c(nrow(bafsegmented), 8)))
  colnames(copyratio_binnedLogR) = c("Chromosome", "Position", "ratioBAF", "ratioBAFphased", "ratioBAF_alt", "ratioBAFphased_alt", "ratioBAFseg", "ratioBAFseg_alt")
  copyratio_binnedLogR$Chromosome = bafsegmented$Chromosome
  copyratio_binnedLogR$Position = bafsegmented$Position

  print("Calculating copy ratios..")
  for (chrom in unique(bafsegmented$Chromosome)) {
    print(chrom)

    baf_chrom = bafsegmented[bafsegmented$Chromosome==chrom,]
    logrseg_chrom = logrsegmented[logrsegmented$Chromosome==chrom,]

    baf_sel = baf_chrom$Position %in% intersect(baf_chrom$Position, logrseg_chrom$Position)
    logrseg_sel = logrseg_chrom$Position %in% intersect(baf_chrom$Position, logrseg_chrom$Position)
    ratio_sel = which(copyratio_binnedLogR$Chromosome==chrom)[baf_sel]

    copyratio_binnedLogR$ratioBAFseg[ratio_sel] = (baf_chrom$BAFseg[baf_sel]*(2^logrseg_chrom$logRseg[logrseg_sel]))
    copyratio_binnedLogR$ratioBAFseg_alt[ratio_sel] = (-(baf_chrom$BAFseg[baf_sel]-1)*(2^logrseg_chrom$logRseg[logrseg_sel]))
  }
  
  background = data.frame(y=seq(0,max.plot.cn,0.5))
  
  print("Plotting..")
  if (platform == "WGS") {
  	sel = seq(1, nrow(allelecounts), 100)
  } else {
	sel = rep(T, nrow(allelecounts))
  }

  plot_title = samplename
  copy_ratio = ggplot(allelecounts[sel,]) +
    geom_hline(data=background, mapping=aes(yintercept=y), colour="black", alpha=0.3) +
    geom_point(mapping=aes(x=Position, y=copy_ratio_binned), alpha=0.5, size=0.9, colour="darkgreen") +
    facet_grid(~Chromosome, scales="free_x", space = "free_x") +
    scale_x_continuous(expand=c(0, 0)) +
    ylim(0,max.plot.cn) + ylab("Copy Ratio") +
    ggtitle(plot_title) +
    theme_bw() + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y = element_text(colour="black",size=18,face="plain"),
                       axis.title.y = element_text(colour="black",size=20,face="plain"),
                       strip.text.x = element_text(colour="black",size=16,face="plain"),
                       plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5))

  if (platform == "WGS") {
	  sel = seq(1, nrow(copyratio_binnedLogR), 100)
  } else {
	  sel = rep(T, nrow(copyratio_binnedLogR))
  }
  as_copy_ratio_seg = ggplot(copyratio_binnedLogR[sel,]) +
    geom_hline(data=background, mapping=aes(yintercept=y), colour="black", alpha=0.3) +
    geom_point(mapping=aes(x=Position, y=ratioBAFseg_alt), alpha=0.5, size=0.9, colour="darkblue") +
    geom_point(mapping=aes(x=Position, y=ratioBAFseg), alpha=0.5, size=0.9, colour="purple") +
    facet_grid(~Chromosome, scales="free_x", space = "free_x") +
    scale_x_continuous(expand=c(0, 0)) +
    ylim(0,max.plot.cn) + ylab("AS Copy Ratio - Segm") +
    theme_bw() + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y = element_text(colour="black",size=18,face="plain"),
                       axis.title.y = element_text(colour="black",size=20,face="plain"),
                       strip.text.x = element_text(colour="black",size=16,face="plain"),
                       plot.title = element_text(colour="black",size=36,face="plain"))
  png(outputfile, width=2000, height=750)
  gridExtra::grid.arrange(gridExtra::arrangeGrob(copy_ratio, as_copy_ratio_seg, ncol=1))
  dev.off()
}

#' Plot relative coverage of tumour and normal
#' 
#' @param samplename Name of the sample for the plot title
#' @param allelecounts Combined allele counts of tumour and normal, read in as a data.frame
#' @param outputfile Full path of file where the figure is to be stored
#' @param max.y The max Y-axis value to be plotted
#' @author sd11
#' @export
coverage_plot = function(samplename, allelecounts, outputfile, max.y=4) {
  
  print("Normalising allele counts..")
  allelecounts$tumour = allelecounts$mutCountT1+allelecounts$mutCountT2
  allelecounts$tumour = allelecounts$tumour / median(allelecounts$tumour, na.rm=T)
  allelecounts$normal = allelecounts$mutCountN1+allelecounts$mutCountN2
  allelecounts$normal = allelecounts$normal / median(allelecounts$normal, na.rm=T)
  
  print("Smoothing data..")
  # res = bin_coverage_tumour(allelecounts, binsize=10000)
  # allelecounts$tumour_binned = res$tumour_binned
  allelecounts$tumour_binned = runmed_data(allelecounts$Chromosome, allelecounts$tumour)
  
  # res = bin_coverage_normal(allelecounts, binsize=10000)
  # allelecounts$normal_binned = res$normal_binned
  # rm(res)
  allelecounts$normal_binned = runmed_data(allelecounts$Chromosome, allelecounts$normal)
  allelecounts$Chromosome = factor(allelecounts$Chromosome, levels=gtools::mixedsort(unique(allelecounts$Chromosome)))
  
  background = data.frame(y=seq(0,2,0.5))
  plot_title = samplename
  p = ggplot(allelecounts[seq(1, nrow(allelecounts), 100),]) +
    geom_hline(data=background, mapping=aes(yintercept=y), colour="black", alpha=0.3) +
    geom_point(mapping=aes(x=Position, y=normal_binned), alpha=0.5, size=0.5, colour="darkgreen") +
    facet_grid(~Chromosome, scales="free_x", space = "free_x") +
    scale_x_continuous(expand=c(0, 0)) +
    ylab("Normal") + scale_y_continuous(breaks=c(0:2), limits=c(0,2)) +
    ggtitle(plot_title) +
    theme_bw() + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y = element_text(colour="black",size=18,face="plain"),
                       axis.title.y = element_text(colour="black",size=20,face="plain"),
                       strip.text.x = element_text(colour="black",size=16,face="plain"),
                       plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5))
  
  background = data.frame(y=seq(0,max.y,0.5))
  p3 = ggplot(allelecounts[seq(1, nrow(allelecounts), 100),]) +
    geom_hline(data=background, mapping=aes(yintercept=y), colour="black", alpha=0.3) +
    geom_point(mapping=aes(x=Position, y=tumour_binned), alpha=0.5, size=0.5, colour="darkgreen") +
    facet_grid(~Chromosome, scales="free_x", space = "free_x") +
    scale_x_continuous(expand=c(0, 0)) +
    ylim(0,max.y) + ylab("Tumour") +
    theme_bw() + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y = element_text(colour="black",size=18,face="plain"),
                       axis.title.y = element_text(colour="black",size=20,face="plain"),
                       strip.text.x = element_text(colour="black",size=16,face="plain"),
                       plot.title = element_text(colour="black",size=36,face="plain"))
  png(outputfile, width=2000, height=750)
  gridExtra::grid.arrange(gridExtra::arrangeGrob(p, p3, ncol=1))
  dev.off()
}



