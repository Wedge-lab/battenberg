#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
#'
create.haplotype.plot = function(chrom.position, points.blue, points.red, x.min, x.max, title, xlab, ylab) {
  par(pch=".", cex=1, cex.main=0.8, cex.axis = 0.6, cex.lab=0.7,yaxp=c(-0.05,1.05,6))
  plot(c(x.min,x.max), c(0,1), type="n",, main=title, xlab=xlab, ylab=ylab)
  points(chrom.position, points.blue, col="blue")
  points(chrom.position, points.red, col="red")
}

#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
#'
create.segmented.plot = function(chrom.position, points.red, points.green, x.min, x.max, title, xlab, ylab) {
  par(mar = c(5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, cex.lab = 2)
  plot(c(x.min,x.max), c(0,1), pch=".", type="n", main=title, xlab=xlab, ylab=ylab)
  points(chrom.position, points.red, pch=".", col="red", cex=2)
  points(chrom.position, points.green, pch=19, cex=0.5, col="green")
}

#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
#'
create.baf.plot = function(chrom.position, points.red.blue, points.darkred, points.darkblue, x.min, x.max, title, xlab, ylab) {
  par(mar = c(5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, cex.lab = 2)
  plot(c(x.min,x.max), c(0,1), pch=".", type = "n", main=title, xlab=xlab, ylab=ylab)
  points(chrom.position, points.red.blue, pch=".", col=ifelse(points.red.blue>0.5, "red", "blue"), cex=2)
  points(chrom.position, points.darkred, pch=19, cex=0.5, col="darkred")
  points(chrom.position, points.darkblue, pch=19, cex=0.5, col="darkblue")
}

#' Function that creates the plots for subclonal copy number
#' Note: This is a plot PER chromosome.
#' @noRd
create.subclonal.cn.plot = function(chrom, chrom.position, LogRposke, LogRchr, BAFchr, BAFsegchr, BAFpvalschr, subcloneres, siglevel, x.min, x.max, title, xlab, ylab.logr, ylab.baf) {
  par(mar=c(2.5,2.5,2.5,0.25), cex=0.4, cex.main=1.5, cex.axis=1, cex.lab=1, mfrow=c(2,1))
  plot(c(x.min, x.max), c(-1,1), pch=".", type="n", main=title, xlab=xlab, ylab=ylab.logr)
  points(LogRposke/1000000, LogRchr, pch=".", col="grey")
  plot(c(x.min, x.max), c(0,1), pch=".", type="n", main=title, xlab=xlab, ylab=ylab.baf)
  points(chrom.position, BAFchr, pch=".", col="grey")
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
#' @noRd
create.bb.plot.average = function(bafsegmented, ploidy, rho, goodnessOfFit, pos_min, pos_max, segment_states_min, segment_states_tot, chr.segs) {
  # Plot main frame and title
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit*100),"%",sep="")
  plot(c(1,nrow(bafsegmented)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
  abline(v=0,lty=1,col="lightgrey")
  # Horizontal lines for y=0 to y=5
  abline(h=c(0:5),lty=1,col="lightgrey")
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
    text(tpos,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)
    abline(v=vpos,lty=1,col="lightgrey")
  }
}

#' Make GW CN plot where subclonal CN is represented by two separate states
#' @noRd
create.bb.plot.subclones = function(bafsegmented, subclones, ploidy, rho, goodnessOfFit, pos_min, pos_max, subcl_min, subcl_max, is_subclonal, is_subclonal_maj, is_subclonal_min, chr.segs) {
	par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
	maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit*100),"%",sep="")
	plot(c(1,nrow(bafsegmented)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
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
		text(tpos,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)
		abline(v=vpos,lty=1,col="lightgrey")
	}
}

#' Code extracted from the first plot in runASCAT.
#' Note: This is a temporary function.
#' @noRd
#'
runascat.plot1 = function(figurename, distancematrix, psis, rhos) {
  require(RColorBrewer)
  png(filename = figurename, width = 1000, height = 1000, res = 1000/7)
  par(mar = c(5,5,0.5,0.5), cex=0.75, cex.lab=2, cex.axis=2)
  hmcol = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  image(log(distancematrix), col = hmcol, axes = F, xlab = "Ploidy", ylab = "Aberrant cell fraction")
  axis(1, at = seq(0, 1, by = 1/5), label = seq(1, 6, by = 1))
  axis(2, at = seq(0, 1/1.05, by = 1/3/1.05), label = seq(0.1, 1, by = 3/10))
  points((psis-1)/5,(rhos-0.1)/0.95,col="green",pch="X", cex = 2)
  dev.off()
}

#' Code extracted from the second plot in runASCAT.
#' Note: This is a temporary function.
#' @noRd
#'
runascat.plot2 = function(figurename, ploidy, rho, goodnessoffit, nonaber, nAf, nBf, chroms, logr, bafsegm, chromnames) {
  png(filename = figurename, width = 2000, height = 500, res = 200)
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessoffit),"%", ifelse(nonaber,", non-aberrant",""),sep="")
    plot(c(1,length(nAf)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(nBf,col="blue",pch = "|")
    points(nAf+nBf,col="purple",pch = "|")
    # don't ask me why, but the "|" sign is not centered, so the lines may need to be shifted..
    abline(v=0,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (i in 1:length(chroms)) {
      chrk = chroms[[i]];
      chrk_hetero = intersect(names(logr)[chrk],names(bafsegm))
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,5,chromnames[i], pos = 1, cex = 2)
      abline(v=vpos,lty=1,col="lightgrey")
    }
    dev.off()
}

#' Code extracted from the second plot in runASCAT.
#' Note: This is a temporary function.
#' @noRd
#'
runascat.plot3 = function(figurename, ploidy, rho, goodnessoffit, nonaber, nA, nB, chroms, logr, bafsegm, chromnames) {
  png(filename = figurename, width = 2000, height = 500, res = 200)
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessoffit),"%", ifelse(nonaber,", non-aberrant",""),sep="")
  plot(c(1,length(nA)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
  points(nA-0.1,col="red",pch = "|")
  points(nB+0.1,col="green",pch = "|")
  # don't ask me why, but the "|" sign is not centered, so the lines may need to be shifted..
  abline(v=0,lty=1,col="lightgrey")
  chrk_tot_len = 0
  for (i in 1:length(ch)) {
    chrk = chroms[[i]];
    chrk_hetero = intersect(names(logr)[chrk],names(bafsegm))
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,5,chomnames[i], pos = 1, cex = 2)
    abline(v=vpos,lty=1,col="lightgrey")
    }
    dev.off()
}

#' Code extracted from the first plot in clonal_ascat runASCAT.
#' Note: This is a temporary function.
#' @noRd
#'
clonal_runascat.plot1 = function(figurename, minim, distmat, psis, rhos) {
  png(filename = figurename, width = 1000, height = 1000, res = 1000/7)
  par(mar = c(5,5,0.5,0.5), cex=0.75, cex.lab=2, cex.axis=2)
  if(minim){ #DCW 240314 reverse colour palette, so blue always corresponds to best region
    hmcol = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  } else {
  hmcol = colorRampPalette(brewer.pal(10, "RdBu"))(256)
  }  
  image(log(distmat), col = hmcol, axes = F, xlab = "Ploidy", ylab = "Aberrant cell fraction")
  axis(1, at = seq(0, 4/4.4, by = 1/4.4), label = seq(1, 5, by = 1))
  axis(2, at = seq(0, 1/1.05, by = 1/3/1.05), label = seq(0.1, 1, by = 3/10))
  points((psis-1)/4.4,(rhos-0.1)/0.95,col="green",pch="X", cex = 2)
  dev.off()
}

#' Code extracted from the second plot in clonal_ascat runASCAT.
#' Note: This is a temporary function and VERY similar to clonal_runascat.plot3
#' @noRd
#'
clonal_runascat.plot2 = function(copynumberprofilespng, rho_opt1, goodnessOfFit_opt1, ploidy_opt1, nA, nB, ch, lrr, bafsegmented, rConf, bConf, confidence) {
  png(filename = copynumberprofilespng, width = 2000, height = 1000, res = 200)
  par(mar = c(0.5,5,5,0.5), mfrow=c(2,1), cex = 0.4, cex.main=3, cex.axis = 2.5)
  maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy_opt1),", aberrant cell fraction: ",sprintf("%2.0f",rho_opt1*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit_opt1),"%",sep="")
  plot(c(1,length(nA)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
  points(nA-0.1,col="red",pch = "|")
  points(nB+0.1,col="green",pch = "|")
  abline(v=0,lty=1,col="lightgrey")
  chrk_tot_len = 0
  for (i in 1:length(ch)) {
    chrk = ch[[i]];
    chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)
    abline(v=vpos,lty=1,col="lightgrey")
  }  
  maintitle = paste("Aberration reliability score (%), average: ", sprintf("%2.0f",mean(confidence,na.rm=T)),"%",sep="")
  plot(c(1,length(nA)), c(0,100), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
  points(confidence,col="blue",pch = "|")
  abline(v=0,lty=1,col="lightgrey")
  chrk_tot_len = 0
  for (i in 1:length(ch)) {
    chrk = ch[[i]];
    chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)
    abline(v=vpos,lty=1,col="lightgrey")
  }
dev.off()
}

#' Code extracted from the third plot in clonal_ascat runASCAT.
#' Note: This is a temporary function and VERY similar to clonal_runascat.plot2
#' @noRd
#'
clonal_runascat.plot3 = function(nonroundedprofilepng, rho_opt1, goodnessOfFit_opt1, ploidy_opt1, nAfull, nBfull, ch, lrr, bafsegmented) {
  png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200)
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy_opt1),", aberrant cell fraction: ",sprintf("%2.0f",rho_opt1*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit_opt1),"%",sep="")
  plot(c(1,length(nAfull)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
  points(nBfull,col="blue",pch = "|")
  points(nAfull+nBfull,col="purple",pch = "|")
  abline(v=0,lty=1,col="lightgrey")
  chrk_tot_len = 0
  for (i in 1:length(ch)) {
    chrk = ch[[i]];
    chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)
    abline(v=vpos,lty=1,col="lightgrey")
  }
  dev.off()
}

