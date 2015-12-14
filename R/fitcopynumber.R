#' Fit copy number
#'
#' Function that will fit a clonal copy number profile to segmented data. It first
#' matches the raw LogR with the segmented BAF to create segmented LogR. Then ASCAT
#' is run to obtain a clonal copy number profile. Beyond logRsegmented it produces 
#' the rho_and_psi file and the cellularity_ploidy file.
#' @param samplename Samplename used to name the segmented logr output file
#' @param outputfile.prefix Prefix used for all output file names, except logRsegmented
#' @param inputfile.baf.segmented Filename that points to the BAF segmented data
#' @param inputfile.baf Filename that points to the raw BAF data
#' @param inputfile.logr Filename that points to the raw LogR data
#' @param dist_choice The distance metric that is used internally to rank clonal copy number solutions
#' @param ascat_dist_choice The distance metric used to obtain an initial cellularity and ploidy estimate
#' @param min.ploidy The minimum ploidy to consider (Default 1.6)
#' @param max.ploidy The maximum ploidy to consider (Default 4.8)
#' @param min.rho The minimum cellularity to consider (Default 0.1)
#' @param min.goodness The minimum goodness of fit for a solution to have to be considered (Default 63)
#' @param uninformative_BAF_threshold The threshold beyond which BAF becomes uninformative (Default 0.51)
#' @param gamma_param Technology parameter, compaction of Log R profiles. Expected decrease in case of deletion in diploid sample, 100 % aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays (Default 1)
#' @param use_preset_rho_psi Boolean whether to use user specified rho and psi values (Default F)
#' @param preset_rho A user specified rho to fit a copy number profile to (Default NA)
#' @param preset_psi A user specified psi to fit a copy number profile to (Default NA)
#' @param read_depth Legacy parameter that is no longer used (Default 30)
#' @author dw9, sd11
#' @export
fit.copy.number = function(samplename, outputfile.prefix, inputfile.baf.segmented, inputfile.baf, inputfile.logr, dist_choice, ascat_dist_choice, min.ploidy=1.6, max.ploidy=4.8, min.rho=0.1, min.goodness=63, uninformative_BAF_threshold=0.51, gamma_param=1, use_preset_rho_psi=F, preset_rho=NA, preset_psi=NA, read_depth=30) {
  
  assert.file.exists(inputfile.baf.segmented)
  assert.file.exists(inputfile.baf)
  assert.file.exists(inputfile.logr)
  
  # Read in the required data
  segmented.BAF.data = read.table(inputfile.baf.segmented, sep="\t", header=T, stringsAsFactors=F) #, row.names=1
  raw.BAF.data = read.table(inputfile.baf, sep="\t", header=T, stringsAsFactors=F)
  raw.logR.data = read.table(inputfile.logr, sep="\t", header=T, stringsAsFactors=F)

  # TODO: strip row names from the input files and check for column names Chromosome, and Position
  
  # Drop NAs
  raw.BAF.data = raw.BAF.data[!is.na(raw.BAF.data[,3]),]
  raw.logR.data = raw.logR.data[!is.na(raw.logR.data[,3]),]
  
  # Chromosome names are sometimes 'chr1', etc.
  if(length(grep("chr",raw.BAF.data[1,1]))>0){
    raw.BAF.data[,1] = gsub("chr","",raw.BAF.data[,1])
  }
  if(length(grep("chr",raw.logR.data[1,1]))>0){
    raw.logR.data[,1] = gsub("chr","",raw.logR.data[,1])
  }
  
  BAF.data = NULL
  logR.data = NULL
  segmented.logR.data = NULL
  matched.segmented.BAF.data = NULL
  chr.names = unique(segmented.BAF.data[,1])
  
  # For each chromosome 
  for(chr in chr.names){
    chr.BAF.data = raw.BAF.data[raw.BAF.data$Chromosome==chr,]
    # Skip the rest if there is no data for this chromosome
    if(nrow(chr.BAF.data)==0){ next }
    # Match segments with chromosome position
    # TODO: Make sure there are column names for this file at the output of segmentation
    chr.segmented.BAF.data = segmented.BAF.data[segmented.BAF.data[,1]==chr,]
    indices = match(chr.segmented.BAF.data[,2],chr.BAF.data$Position )

    if (sum(is.na(indices))==length(indices) | length(indices)==0) {
	    next
    }
    
    # Drop NAs here too
    chr.segmented.BAF.data = chr.segmented.BAF.data[!is.na(indices),]
    
    # Append the segmented data
    matched.segmented.BAF.data = rbind(matched.segmented.BAF.data, chr.segmented.BAF.data)
    BAF.data = rbind(BAF.data, chr.BAF.data[indices[!is.na(indices)],])
    
    # Append raw LogR
    chr.logR.data = raw.logR.data[raw.logR.data$Chromosome==chr,]
    indices = match(chr.segmented.BAF.data[,2],chr.logR.data$Position)
    logR.data = rbind(logR.data, chr.logR.data[indices[!is.na(indices)],])
    chr.segmented.logR.data = chr.logR.data[indices[!is.na(indices)],]
    
    # Append segmented LogR
    segs = rle(chr.segmented.BAF.data[,5])$lengths
    cum.segs = c(0,cumsum(segs))
    for(s in 1:length(segs)){
      chr.segmented.logR.data[(cum.segs[s]+1):cum.segs[s+1],3] = mean(chr.segmented.logR.data[(cum.segs[s]+1):cum.segs[s+1],3])
    }
    segmented.logR.data = rbind(segmented.logR.data,chr.segmented.logR.data)
  }
  names(matched.segmented.BAF.data)[5] = samplename
  row.names(segmented.logR.data) = row.names(matched.segmented.BAF.data)
  row.names(logR.data) = row.names(matched.segmented.BAF.data)
  
  write.table(segmented.logR.data,paste(samplename,".logRsegmented.txt",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
  
  segBAF = 1-matched.segmented.BAF.data[,5]
  segLogR = segmented.logR.data[,3]
  logR = logR.data[,3]
  names(segBAF) = rownames(matched.segmented.BAF.data)
  names(segLogR) = rownames(matched.segmented.BAF.data)
  names(logR) = rownames(matched.segmented.BAF.data)
  print(unique(logR.data[,1]))
  chr.segs = NULL
  for(ch in 1:length(chr.names)){
    chr.segs[[ch]] = which(logR.data[,1]==chr.names[ch])
  }
  
  if(use_preset_rho_psi){
    ascat_optimum_pair = list(rho=preset_rho, psi = preset_psi, ploidy = preset_psi)
  }else{
    distance.outfile=paste(outputfile.prefix,"distance.png",sep="",collapse="") # kjd 20-2-2014
    copynumberprofile.outfile=paste(outputfile.prefix,"copynumberprofile.png",sep="",collapse="") # kjd 20-2-2014
    nonroundedprofile.outfile=paste(outputfile.prefix,"nonroundedprofile.png",sep="",collapse="") # kjd 20-2-2014
    
    ascat_optimum_pair=runASCAT(logR, 1-BAF.data[,3], segLogR, segBAF, chr.segs, ascat_dist_choice,distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, gamma=gamma_param, allow100percent=T, reliabilityFile=NA, min.ploidy, max.ploidy, min.rho, min.goodness) # kjd 4-2-2014
  }
  
  distance.outfile=paste(outputfile.prefix,"second_distance.png",sep="",collapse="") # kjd 20-2-2014
  copynumberprofile.outfile=paste(outputfile.prefix,"second_copynumberprofile.png",sep="",collapse="") # kjd 20-2-2014
  nonroundedprofile.outfile=paste(outputfile.prefix,"second_nonroundedprofile.png",sep="",collapse="") # kjd 20-2-2014
  
  # All is set up, now run ASCAT to obtain a clonal copynumber profile
  out = run_clonal_ASCAT( logR, 1-BAF.data[,3], segLogR, segBAF, chr.segs, matched.segmented.BAF.data, ascat_optimum_pair, dist_choice, distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, gamma_param=gamma_param, read_depth, uninformative_BAF_threshold, allow100percent=T, reliabilityFile=NA) # kjd 21-2-2014
  
  ascat_optimum_pair_fraction_of_genome = out$output_optimum_pair_without_ref
  ascat_optimum_pair_ref_seg = out$output_optimum_pair
  is.ref.better = out$is.ref.better

  # Save rho, psi and ploidy for future reference
  rho_psi_output = data.frame(rho = c(ascat_optimum_pair$rho,ascat_optimum_pair_fraction_of_genome$rho,ascat_optimum_pair_ref_seg$rho),psi = c(ascat_optimum_pair$psi,ascat_optimum_pair_fraction_of_genome$psi,ascat_optimum_pair_ref_seg$psi), ploidy = c(ascat_optimum_pair$ploidy,ascat_optimum_pair_fraction_of_genome$ploidy,ascat_optimum_pair_ref_seg$ploidy), distance = c(NA,out$distance_without_ref,out$distance), is.best = c(NA,!is.ref.better,is.ref.better),row.names=c("ASCAT","FRAC_GENOME","REF_SEG"))
  write.table(rho_psi_output,paste(outputfile.prefix,"rho_and_psi.txt",sep=""),quote=F,sep="\t")

  # Create user friendly cellularity and ploidy output file
  cellularity_ploidy_output = data.frame(cellularity = c(ascat_optimum_pair_fraction_of_genome$rho), ploidy = c(ascat_optimum_pair_fraction_of_genome$ploidy), psi = c(ascat_optimum_pair_fraction_of_genome$psi))
  write.table(cellularity_ploidy_output, paste(outputfile.prefix,"cellularity_ploidy.txt",sep=""), quote=F, sep="\t", row.names=F)
}

#' Fit subclonal copy number
#' 
#' This function fits a subclonal copy number profile where a clonal profile is unlikely.
#' It goes over each segment of a clonal copy number profile and does a simple t-test. If the
#' test is significant it is unlikely that the data can be explained by a single copy number 
#' state. We therefore fit a second state, i.e. there are two cellular populations with each
#' a different state: Subclonal copy number.
#' @param sample.name Name of the sample, used in figures
#' @param baf.segmented.file String that points to a file with segmented BAF output
#' @param logr.file String that points to the raw LogR file to be used in the subclonal copy number figures
#' @param rho.psi.file String pointing to the rho_and_psi file generated by fit.copy.number
#' @param output.file Filename of the file where the final copy number fit will be written to
#' @param output.figures.prefix Prefix of the filenames for the chromosome specific copy number figures
#' @param output.gw.figures.prefix Prefix of the filenames for the genome wide copy number figures
#' @param chr_names Vector of allowed chromosome names
#' @param gamma Technology specific scaling parameter for LogR (Default 1)
#' @param segmentation.gamma Legacy parameter that is no longer used (Default NA)
#' @param siglevel Threshold under which a p-value becomes significant. When it is significant a second copy number state will be fitted (Default 0.05)
#' @param maxdist Slack in BAF space to allow a segment to be off it's optimum before becoming significant. A segment becomes significant very quickly when a breakpoint is missed, this parameter alleviates the effect (Default 0.01)
#' @param noperms The number of permutations to be run when bootstrapping the confidence intervals on the copy number state of each segment (Default 1000)
#' @author dw9, sd11
#' @export
callSubclones = function(sample.name, baf.segmented.file, logr.file, rho.psi.file, output.file, output.figures.prefix, output.gw.figures.prefix, chr_names, gamma=1, segmentation.gamma=NA, siglevel=0.05, maxdist=0.01, noperms=1000) {
  
  # Load rho/psi/goodness of fit
  res = load.rho.psi.file(rho.psi.file)
  rho = res$rho
  psit = res$psit
  psi = rho*psit + 2 * (1-rho) # psi of all cells
  goodness = res$goodness
  
  # Load the BAF segmented data
  BAFvals = read.table(baf.segmented.file, sep="\t", header=T, stringsAsFactors=F) #, row.names=F
  if (colnames(BAFvals)[1] == "X") {
	  # If there were rownames, then delete this column. Should not be an issue with new BB runs
	  BAFvals = BAFvals[,-1]
  }

  BAF = BAFvals[,3]
  BAFphased = BAFvals[,4]
  BAFseg = BAFvals[,5]
  
  # Save SNP positions separately
  SNPpos = BAFvals[,c(1,2)]
  
  # Load the raw LogR data
  LogRvals = read.table(logr.file,sep="\t", header=T, stringsAsFactors=F)
  if (colnames(LogRvals)[1] == "X") {
	  # If there were rownames, then delete this column. Should not be an issue with new BB runs
	  LogRvals = LogRvals[,-1]
  }

  # Chromosome names are sometimes 'chr1', etc.
  if(length(grep("chr",LogRvals[1,1]))>0){
	      LogRvals[,1] = gsub("chr","",LogRvals[,1])
  }
  
  ctrans = c(1:length(chr_names))
  names(ctrans) = chr_names
  ctrans.logR = c(1:length(chr_names))
  names(ctrans.logR) = chr_names
  
  LogRpos = as.vector(ctrans.logR[as.vector(LogRvals[,1])]*1000000000+LogRvals[,2])
  BAFpos = as.vector(ctrans[as.vector(BAFvals[,1])]*1000000000+BAFvals[,2])
  
  #DCW 240314
  switchpoints = c(0,which(BAFseg[-1] != BAFseg[-(length(BAFseg))] | BAFvals[-1,1] != BAFvals[-nrow(BAFvals),1]),length(BAFseg))
  BAFlevels = BAFseg[switchpoints[-1]]
  
  pval = NULL
  BAFpvals = vector(length=length(BAFseg))
  subcloneres = NULL

  print(head(BAFvals))
  print(head(LogRvals))

  for (i in 1:length(BAFlevels)) {
    l = BAFlevels[i]
    
    # Make sure that BAF>=0.5, otherwise nMajor and nMinor may be the wrong way around
    l = max(l,1-l)
    
    BAFke = BAFphased[(switchpoints[i]+1):switchpoints[i+1]]
    
    #startpos = min(BAFpos[names(BAFke)])
    #endpos = max(BAFpos[names(BAFke)])
    startpos = min(BAFpos[(switchpoints[i]+1):switchpoints[i+1]])
    endpos = max(BAFpos[(switchpoints[i]+1):switchpoints[i+1]])
    chrom = names(ctrans[floor(startpos/1000000000)])
    LogR = mean(LogRvals[LogRpos>=startpos&LogRpos<=endpos & !is.infinite(LogRvals[,3]),3],na.rm=T)
    
    # if we don't have a value for LogR, fill in 0
    if (is.na(LogR)) {
      LogR = 0
    }
    nMajor = (rho-1+l*psi*2^(LogR/gamma))/rho
    nMinor = (rho-1+(1-l)*psi*2^(LogR/gamma))/rho
    
    # Increase nMajor and nMinor together, to avoid impossible combinations (with negative subclonal fractions)
    if (nMinor<0) {
      if (l==1) {
        # Avoid calling infinite copy number
        nMajor = 1000
      } else {
        nMajor = nMajor + l * (0.01 - nMinor) / (1-l)
      }
      nMinor = 0.01  	
    }
    
    # Note that these are sorted in the order of ascending BAF:
    nMaj = c(floor(nMajor), ceiling(nMajor), floor(nMajor), ceiling(nMajor))
    nMin = c(ceiling(nMinor), ceiling(nMinor), floor(nMinor), floor(nMinor))
    x = floor(nMinor)
    y = floor(nMajor)
    
    # Total copy number, to determine priority options
    ntot = nMajor + nMinor
    
    levels = (1-rho+rho*nMaj)/(2-2*rho+rho*(nMaj+nMin))
    # Problem if rho=1 and nMaj=0 and nMin=0
    levels[nMaj==0 & nMin==0] = 0.5
    
    #DCW - just test corners on the nearest edge to determine clonality
    # If the segment is called as subclonal, this is the edge that will be used to determine the subclonal proportions that are reported first
    all.edges = orderEdges(levels, l, ntot,x,y)
    nMaj.test = all.edges[1,c(1,3)]
    nMin.test = all.edges[1,c(2,4)]
    test.levels = (1-rho+rho*nMaj.test)/(2-2*rho+rho*(nMaj.test+nMin.test))
    whichclosestlevel.test = which.min(abs(test.levels-l))
    
    # Test whether a segment should be subclonal
    if (sd(BAFke)==0) {
      pval[i] = 0 # problem caused by segments with constant BAF (usually 1 or 2)
    } else {
      pval[i] = t.test(BAFke, alternative="two.sided", mu=test.levels[whichclosestlevel.test])$p.value
    }
    if (abs(l-test.levels[whichclosestlevel.test])<maxdist) {  
      pval[i] = 1
    }
    
    #DCW 240314
    BAFpvals[(switchpoints[i]+1):switchpoints[i+1]] = pval[i]
    
    # If the difference is significant, call subclonal level
    if (pval[i] <= siglevel) {
      
      all.edges = orderEdges(levels, l, ntot,x,y)
      # Switch order, so that negative copy numbers are at the end
      na.indices = which(is.na(rowSums(all.edges)))
      if (length(na.indices)>0) {
        all.edges = rbind(all.edges[-na.indices,], all.edges[na.indices,])
      }
      nMaj1 = all.edges[,1]
      nMin1 = all.edges[,2]
      nMaj2 = all.edges[,3]
      nMin2 = all.edges[,4]
      
      tau = (1 - rho + rho * nMaj2 - 2 * l * (1 - rho) - l * rho * (nMin2 + nMaj2)) / (l * rho * (nMin1 + nMaj1) - l * rho * (nMin2 + nMaj2) - rho * nMaj1 + rho * nMaj2)
      sdl = sd(BAFke,na.rm=T)/sqrt(sum(!is.na(BAFke)))
      sdtau = abs((1 - rho + rho * nMaj2 - 2 * (l+sdl) * (1 - rho) - (l+sdl) * rho * (nMin2 + nMaj2)) / ((l+sdl) * rho * (nMin1 + nMaj1) - (l+sdl) * rho * (nMin2 + nMaj2) - rho * nMaj1 + rho * nMaj2) - tau) / 2 +
        abs((1 - rho + rho * nMaj2 - 2 * (l-sdl) * (1 - rho) - (l-sdl) * rho * (nMin2 + nMaj2)) / ((l-sdl) * rho * (nMin1 + nMaj1) - (l-sdl) * rho * (nMin2 + nMaj2) - rho * nMaj1 + rho * nMaj2) - tau) / 2

      # Bootstrapping to obtain 95% confidence intervals
      sdtaubootstrap = vector(length=length(tau), mode="numeric")
      tau25 = vector(length=length(tau), mode="numeric")
      tau975 = vector(length=length(tau), mode="numeric")
      
      for (option in 1:length(tau)) {
        nMaj1o = nMaj1[option]
        nMin1o = nMin1[option]
        nMaj2o = nMaj2[option]
        nMin2o = nMin2[option]
        set.seed(as.integer(Sys.time()))
        permFraction = vector(length=noperms,mode="numeric")
        for (j in 1:noperms) {
          permBAFs=sample(BAFke,length(BAFke),replace=T)
          permMeanBAF=mean(permBAFs)
          permFraction[j] = (1 - rho + rho * nMaj2o - 2 * permMeanBAF * (1 - rho) - permMeanBAF * rho * (nMin2o + nMaj2o)) / (permMeanBAF * rho * (nMin1o + nMaj1o) - permMeanBAF * rho * (nMin2o + nMaj2o) - rho * nMaj1o + rho * nMaj2o)
        }
        orderedFractions = sort(permFraction)
        sdtaubootstrap[option] = sd(permFraction)
        tau25[option] = orderedFractions[25]
        tau975[option] = orderedFractions[975]
      }
      
      subcloneres = rbind(subcloneres,c(chrom,startpos-floor(startpos/1000000000)*1000000000,
                                        endpos-floor(endpos/1000000000)*1000000000,l,pval[i],LogR,ntot,
                                        nMaj1[1],nMin1[1],tau[1],nMaj2[1],nMin2[1],1-tau[1],sdtau[1],sdtaubootstrap[1],tau25[1],tau975[1],
                                        nMaj1[2],nMin1[2],tau[2],nMaj2[2],nMin2[2],1-tau[2],sdtau[2],sdtaubootstrap[2],tau25[2],tau975[2],
                                        nMaj1[3],nMin1[3],tau[3],nMaj2[3],nMin2[3],1-tau[3],sdtau[3],sdtaubootstrap[3],tau25[3],tau975[3],
                                        nMaj1[4],nMin1[4],tau[4],nMaj2[4],nMin2[4],1-tau[4],sdtau[4],sdtaubootstrap[4],tau25[4],tau975[4],
                                        nMaj1[5],nMin1[5],tau[5],nMaj2[5],nMin2[5],1-tau[5],sdtau[5],sdtaubootstrap[5],tau25[5],tau975[5],
                                        nMaj1[6],nMin1[6],tau[6],nMaj2[6],nMin2[6],1-tau[6],sdtau[6],sdtaubootstrap[6],tau25[6],tau975[6]))
    }else {
      #if called as clonal, use the best corner from the nearest edge 
      subcloneres = rbind(subcloneres,c(chrom,startpos-floor(startpos/1000000000)*1000000000,
                                        endpos-floor(endpos/1000000000)*1000000000,l,pval[i],LogR,ntot,
                                        nMaj.test[whichclosestlevel.test],nMin.test[whichclosestlevel.test],1,rep(NA,57)))
      
    }
  }

  colnames(subcloneres) = c("chr","startpos","endpos","BAF","pval","LogR","ntot",
                            "nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","frac2_A","SDfrac_A","SDfrac_A_BS","frac1_A_0.025","frac1_A_0.975",
                            "nMaj1_B","nMin1_B","frac1_B","nMaj2_B","nMin2_B","frac2_B","SDfrac_B","SDfrac_B_BS","frac1_B_0.025","frac1_B_0.975",
                            "nMaj1_C","nMin1_C","frac1_C","nMaj2_C","nMin2_C","frac2_C","SDfrac_C","SDfrac_C_BS","frac1_C_0.025","frac1_C_0.975",
                            "nMaj1_D","nMin1_D","frac1_D","nMaj2_D","nMin2_D","frac2_D","SDfrac_D","SDfrac_D_BS","frac1_D_0.025","frac1_D_0.975",
                            "nMaj1_E","nMin1_E","frac1_E","nMaj2_E","nMin2_E","frac2_E","SDfrac_E","SDfrac_E_BS","frac1_E_0.025","frac1_E_0.975",
                            "nMaj1_F","nMin1_F","frac1_F","nMaj2_F","nMin2_F","frac2_F","SDfrac_F","SDfrac_F_BS","frac1_F_0.025","frac1_F_0.975")
  
  write.table(subcloneres, output.file, quote=F, col.names=T, row.names=F, sep="\t")

  # Create a plot per chromosome that shows the segments with their CN state in text
  for (chr in chr_names) {
    pos = SNPpos[SNPpos[,1]==chr, 2]
    #if no points to plot, skip
    if (length(pos)==0) { next }

    png(filename = paste(output.figures.prefix, chr,".png",sep=""), width = 2000, height = 2000, res = 200)
    create.subclonal.cn.plot(chrom=chr,
                             chrom.position=pos/1000000, 
                             LogRposke=LogRvals[LogRvals[,1]==chr,2], 
                             LogRchr=LogRvals[LogRvals[,1]==chr,3], 
                             BAFchr=BAF[SNPpos[,1]==chr], 
                             BAFsegchr=BAFseg[SNPpos[,1]==chr], 
                             BAFpvalschr=BAFpvals[SNPpos[,1]==chr],
                             subcloneres=subcloneres, 
                             siglevel=siglevel, 
                             x.min=min(pos)/1000000, 
                             x.max=max(pos)/1000000, 
                             title=paste(sample.name,", chromosome ", chr, sep=""), 
                             xlab="Position (Mb)", 
                             ylab.logr="LogR", 
                             ylab.baf="BAF (phased)")
    dev.off()
  }
 
  # Cast columns back to numeric
  subclones = as.data.frame(subcloneres)
  subclones[,2:ncol(subclones)] = sapply(2:ncol(subclones), function(x) { as.numeric(as.character(subclones[,x])) })
  plot.gw.subclonal.cn(subclones=subclones, BAFvals=BAFvals, rho=rho, psi=psi, goodness=goodness, output.gw.figures.prefix=output.gw.figures.prefix, chr.names=chr_names)

}

#' Plot the copy number genome wide in two different ways. This creates the Battenberg average 
#' profile where subclonal copy number is represented as a mixture of two different states and 
#' the Battenberg subclones profile where subclonal copy number is plotted as two different
#' separate states. The thickness of the line represents the fraction of tumour cells carying
#' the particular state.
#' @noRd
plot.gw.subclonal.cn = function(subclones, BAFvals, rho, psi, goodness, output.gw.figures.prefix, chr.names) {
  # Map start and end of each segment into the BAF values. The plot uses the index of this BAF table as x-axis
  pos_min = array(NA, nrow(subclones))
  pos_max = array(NA, nrow(subclones))
  for (i in 1:nrow(subclones)) {
    segm_chr = subclones$chr[i] == BAFvals$Chromosome & subclones$startpos[i] < BAFvals$Position & subclones$endpos[i] >= BAFvals$Position
    pos_min[i] = min(which(segm_chr))
    pos_max[i] = max(which(segm_chr))
  }
  
  # For those segments that are subclonal, Obtain the second state.
  is_subclonal = which(subclones$frac1_A < 1)
  subcl_min = array(NA, length(is_subclonal))
  subcl_max = array(NA, length(is_subclonal))
  for (i in 1:length(is_subclonal)) {
    segment_index = is_subclonal[i]
    segm_chr = subclones$chr[segment_index] == BAFvals$Chromosome & subclones$startpos[segment_index] < BAFvals$Position & subclones$endpos[segment_index] >= BAFvals$Position
    subcl_min[i] = min(which(segm_chr))
    subcl_max[i] = max(which(segm_chr))
  }
  
  # Determine whether it's the major or the minor allele that is represented by two states
  is_subclonal_maj = abs(subclones$nMaj1_A - subclones$nMaj2_A) > 0
  is_subclonal_min = abs(subclones$nMin1_A - subclones$nMin2_A) > 0
  is_subclonal_maj[is.na(is_subclonal_maj)] = F
  is_subclonal_min[is.na(is_subclonal_min)] = F
  
  # BB represents subclonal CN as a mixture of two CN states. Calculate this mixture for both minor allele and total CN.
  #segment_states_min = subclones$nMin1_A * ifelse(is_subclonal_min, subclones$frac1_A, 1)  + ifelse(is_subclonal_min, subclones$nMin2_A, 0) * ifelse(is_subclonal_min, subclones$frac2_A, 0) 
  #segment_states_tot = (subclones$nMaj1_A+subclones$nMin1_A) * ifelse(is_subclonal_maj, subclones$frac1_A, 1) + ifelse(is_subclonal_maj, subclones$nMaj2_A+subclones$nMin2_A, 0) * ifelse(is_subclonal_maj, subclones$frac2_A, 0) 

  segment_states_min = subclones$nMin1_A * ifelse(is_subclonal_min, subclones$frac1_A, 1)  + ifelse(is_subclonal_min, subclones$nMin2_A, 0) * ifelse(is_subclonal_min, subclones$frac2_A, 0) 
  segment_states_maj = subclones$nMaj1_A * ifelse(is_subclonal_maj, subclones$frac1_A, 1)  + ifelse(is_subclonal_maj, subclones$nMaj2_A, 0) * ifelse(is_subclonal_maj, subclones$frac2_A, 0) 
  segment_states_tot = segment_states_maj + segment_states_min
  
  # Determine which SNPs are on which chromosome, to be used as a proxy for chromosome size in the plots
  chr.segs = lapply(1:length(chr.names), function(ch) { which(BAFvals$Chromosome==chr.names[ch]) })
  
  # Plot subclonal copy number as mixtures of two states
  png(filename = paste(output.gw.figures.prefix, "_average.png", sep=""), width = 2000, height = 500, res = 200)
  create.bb.plot.average(bafsegmented=BAFvals, 
                            ploidy=psi, 
                            rho=rho, 
                            goodnessOfFit=goodness, 
                            pos_min=pos_min, 
                            pos_max=pos_max, 
                            segment_states_min=segment_states_min, 
                            segment_states_tot=segment_states_tot, 
                            chr.segs=chr.segs)
  dev.off()
  
  # Plot subclonal copy number as two separate states
  png(filename = paste(output.gw.figures.prefix, "_subclones.png", sep=""), width = 2000, height = 500, res = 200)
  create.bb.plot.subclones(bafsegmented=BAFvals, 
                             subclones=subclones,
                             ploidy=psi, 
                             rho=rho, 
                             goodnessOfFit=goodness, 
                             pos_min=pos_min, 
                             pos_max=pos_max, 
                             subcl_min=subcl_min, 
                             subcl_max=subcl_max, 
                             is_subclonal=is_subclonal,
			     is_subclonal_maj=is_subclonal_maj,
			     is_subclonal_min=is_subclonal_min,
			     chr.segs=chr.segs)
  dev.off()
}

#' Load the rho and psi estimates from a file.
#' @noRd
load.rho.psi.file = function(rho.psi.file) {
  rho_psi_info = read.table(rho.psi.file, header=T, sep="\t", stringsAsFactors=F)
  # Always use best solution from grid search - reference segment sometimes gives strange results
  rho = rho_psi_info$rho[rownames(rho_psi_info)=="FRAC_GENOME"] # rho = tumour percentage (called tp in previous versions)
  psit = rho_psi_info$psi[rownames(rho_psi_info)=="FRAC_GENOME"] # psi of tumour cells
  goodness = rho_psi_info$distance[rownames(rho_psi_info)=="FRAC_GENOME"] # goodness of fit
  return(list(rho=rho, psit=psit, goodness=goodness))
}
