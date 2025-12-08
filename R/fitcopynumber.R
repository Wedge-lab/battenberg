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
#' @param max.rho The maximum cellularity to consider (Default 1.0)
#' @param min.goodness The minimum goodness of fit for a solution to have to be considered (Default 63)
#' @param uninformative_BAF_threshold The threshold beyond which BAF becomes uninformative (Default 0.51)
#' @param gamma_param Technology parameter, compaction of Log R profiles. Expected decrease in case of deletion in diploid sample, 100 "\%" aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays (Default 1)
#' @param use_preset_rho_psi Boolean whether to use user specified rho and psi values (Default F)
#' @param preset_rho A user specified rho to fit a copy number profile to (Default NA)
#' @param preset_psi A user specified psi to fit a copy number profile to (Default NA)
#' @param read_depth Legacy parameter that is no longer used (Default 30)
#' @param analysis A String representing the type of analysis to be run, this determines whether the distance figure is produced (Default paired)
#' @author dw9, sd11
#' @export
fit.copy.number = function(samplename, outputfile.prefix, inputfile.baf.segmented, inputfile.baf, inputfile.logr, dist_choice, ascat_dist_choice, min.ploidy=1.6, max.ploidy=4.8, min.rho=0.1,  max.rho=1.0, min.goodness=63, uninformative_BAF_threshold=0.51, gamma_param=1, use_preset_rho_psi=F, preset_rho=NA, preset_psi=NA, read_depth=30, analysis="paired", nthreads, enhanced_grid_search=F) {
  
  assert.file.exists(inputfile.baf.segmented)
  assert.file.exists(inputfile.baf)
  assert.file.exists(inputfile.logr)
  # Check for enough options supplied for rho and psi
  if ((max.ploidy - min.ploidy) < 0.05) {
    stop(paste("Supplied ploidy range must be larger than 0.05: ", min.ploidy, "-", max.ploidy, sep=""))
  }
  if ((max.rho - min.rho) < 0.01) {
    stop(paste("Supplied rho range must be larger than 0.01: ", min.rho, "-", max.rho, sep=""))
  }
  
  # Read in the required data
  segmented.BAF.data = as.data.frame(read_bafsegmented(inputfile.baf.segmented))
  raw.BAF.data = as.data.frame(read_baf(inputfile.baf))
  raw.logR.data = as.data.frame(read_logr(inputfile.logr))
  
  # Assign rownames as those are required by various clonal_ascat.R functions
  # If there are duplicates (possible with old versions of BB) then remove those
  identifiers = paste(segmented.BAF.data[,1], segmented.BAF.data[,2], sep="_")
  dups = which(duplicated(identifiers))
  if (length(dups) > 0) {
    segmented.BAF.data = segmented.BAF.data[-dups,]
    identifiers = identifiers[-dups]
  }
  rownames(segmented.BAF.data) = identifiers
  
  # Drop NAs
  raw.BAF.data = raw.BAF.data[!is.na(raw.BAF.data[,3]),]
  raw.logR.data = raw.logR.data[!is.na(raw.logR.data[,3]),]
  
  ## Chromosome names are sometimes 'chr1', etc.
  #if(length(grep("chr",raw.BAF.data[1,1]))>0){
  #  raw.BAF.data[,1] = gsub("chr","",raw.BAF.data[,1])
  #}
  #if(length(grep("chr",raw.logR.data[1,1]))>0){
  #  raw.logR.data[,1] = gsub("chr","",raw.logR.data[,1])
  #}
  
  BAF.data = list()
  logR.data = list()
  segmented.logR.data = list()
  matched.segmented.BAF.data = list()
  gsubchr = function(chr) gsub("chr","",as.character(chr))
  
  chr.names = gsubchr(unique(segmented.BAF.data[,1]))
  
  segmented.BAF.data$Chromosome = gsubchr(segmented.BAF.data$Chromosome)
  raw.BAF.data$Chromosome = gsubchr(raw.BAF.data$Chromosome)
  raw.logR.data$Chromosome =gsubchr(raw.logR.data$Chromosome)
  
  baf_segmented_split = split(segmented.BAF.data, f=segmented.BAF.data$Chromosome)
  baf_split = split(raw.BAF.data, f=raw.BAF.data$Chromosome)
  logr_split = split(raw.logR.data, f=raw.logR.data$Chromosome)
  
  # For each chromosome
  for(chr in chr.names){
    chr.BAF.data = baf_split[[chr]]
    
    # Skip the rest if there is no data for this chromosome
    if(is.null(chr.BAF.data) || nrow(chr.BAF.data)==0){ next }
    # Match segments with chromosome position
    chr.segmented.BAF.data = baf_segmented_split[[chr]]
    indices = match(chr.segmented.BAF.data[,2],chr.BAF.data$Position )
    
    if (sum(is.na(indices))==length(indices) | length(indices)==0) {
      next
    }
    
    # Drop NAs here too
    chr.segmented.BAF.data = chr.segmented.BAF.data[!is.na(indices),]
    
    # Append the segmented data
    matched.segmented.BAF.data[[chr]] = chr.segmented.BAF.data
    BAF.data[[chr]] = chr.BAF.data[indices[!is.na(indices)],]
    
    # Append raw LogR
    chr.logR.data = logr_split[[chr]]
    indices = match(chr.segmented.BAF.data[,2],chr.logR.data$Position)
    logR.data[[chr]] = chr.logR.data[indices[!is.na(indices)],]
    chr.segmented.logR.data = chr.logR.data[indices[!is.na(indices)],]
    
    # Append segmented LogR
    segs = rle(chr.segmented.BAF.data[,5])$lengths
    cum.segs = c(0,cumsum(segs))
    for(s in 1:length(segs)){
      chr.segmented.logR.data[(cum.segs[s]+1):cum.segs[s+1],3] = mean(chr.segmented.logR.data[(cum.segs[s]+1):cum.segs[s+1],3], na.rm=T)
    }
    segmented.logR.data[[chr]] = chr.segmented.logR.data
  }
  
  # Sync the dataframes
  selection = c()
  for (chrom in chr.names) {
    matched.segmented.BAF.data.chr = matched.segmented.BAF.data[[chrom]] #matched.segmented.BAF.data[matched.segmented.BAF.data[,1]==chrom,]
    logR.data.chr = logR.data[[chrom]] #logR.data[logR.data[,1]==chrom,]
    
    selection = matched.segmented.BAF.data.chr[,2] %in% logR.data.chr[,2]
    matched.segmented.BAF.data[[chrom]] = matched.segmented.BAF.data.chr[selection,]
    segmented.logR.data[[chrom]] = segmented.logR.data[[chrom]][selection,]
  }
  
  # Combine the split data frames into a single for the subsequent steps
  matched.segmented.BAF.data = do.call(rbind, matched.segmented.BAF.data)
  segmented.logR.data = do.call(rbind, segmented.logR.data)
  BAF.data = do.call(rbind, BAF.data)
  logR.data = do.call(rbind, logR.data)
  names(matched.segmented.BAF.data)[5] = samplename
  
  # write out the segmented logR data
  row.names(segmented.logR.data) = row.names(matched.segmented.BAF.data)
  row.names(logR.data) = row.names(matched.segmented.BAF.data)
  write.table(segmented.logR.data,paste(samplename,".logRsegmented.txt",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
  
  # Prepare the data for going into the runASCAT functions
  segBAF = 1-matched.segmented.BAF.data[,5]
  segLogR = segmented.logR.data[,3]
  logR = logR.data[,3]
  names(segBAF) = rownames(matched.segmented.BAF.data)
  names(segLogR) = rownames(matched.segmented.BAF.data)
  names(logR) = rownames(matched.segmented.BAF.data)
  
  chr.segs = NULL
  for(ch in 1:length(chr.names)){
    chr.segs[[ch]] = which(logR.data[,1]==chr.names[ch])
  }
  
  if(use_preset_rho_psi){
    ascat_optimum_pair = list(rho=preset_rho, psi = preset_psi, ploidy = preset_psi)
  }else{
    distance.outfile=paste(outputfile.prefix, "distance.png", sep="", collapse="") # kjd 20-2-2014
    copynumberprofile.outfile=paste(outputfile.prefix, "copynumberprofile.png", sep="", collapse="") # kjd 20-2-2014
    nonroundedprofile.outfile=paste(outputfile.prefix, "nonroundedprofile.png", sep="", collapse="") # kjd 20-2-2014
    cnaStatusFile = paste(outputfile.prefix, "copynumber_solution_status.txt", sep="", collapse="")
    
    if(enhanced_grid_search) {
      ascat_optimum_pair = runASCAT_enhanced(logR, 1-BAF.data[,3], segLogR, segBAF, chr.segs, ascat_dist_choice,distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, cnaStatusFile=cnaStatusFile, gamma=gamma_param, allow100percent=T, reliabilityFile=NA, min.ploidy=min.ploidy, max.ploidy=max.ploidy, min.rho=min.rho, max.rho=max.rho, min.goodness, chr.names=chr.names, analysis=analysis, uninformative_BAF_threshold=uninformative_BAF_threshold, verbose=TRUE)
    } else {
      ascat_optimum_pair = runASCAT(logR, 1-BAF.data[,3], segLogR, segBAF, chr.segs, ascat_dist_choice,distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, cnaStatusFile=cnaStatusFile, gamma=gamma_param, allow100percent=T, reliabilityFile=NA, min.ploidy=min.ploidy, max.ploidy=max.ploidy, min.rho=min.rho, max.rho=max.rho, min.goodness, chr.names=chr.names, analysis=analysis) # kjd 4-2-2014
    }
  }
  
  distance.outfile=paste(outputfile.prefix,"second_distance.png",sep="",collapse="") # kjd 20-2-2014
  copynumberprofile.outfile=paste(outputfile.prefix,"second_copynumberprofile.png",sep="",collapse="") # kjd 20-2-2014
  nonroundedprofile.outfile=paste(outputfile.prefix,"second_nonroundedprofile.png",sep="",collapse="") # kjd 20-2-2014
  
  # All is set up, now run ASCAT to obtain a clonal copynumber profile
  out = run_clonal_ASCAT( logR, 1-BAF.data[,3], segLogR, segBAF, chr.segs, matched.segmented.BAF.data, ascat_optimum_pair, dist_choice, distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, gamma_param=gamma_param, read_depth, uninformative_BAF_threshold, allow100percent=T, reliabilityFile=NA,  psi_min_initial=min.ploidy, psi_max_initial=max.ploidy, rho_min_initial=min.rho, rho_max_initial=max.rho, chr.names=chr.names) # kjd 21-2-2014
  
  ascat_optimum_pair_fraction_of_genome = out$output_optimum_pair_without_ref
  ascat_optimum_pair_ref_seg = out$output_optimum_pair
  is.ref.better = out$is.ref.better
  
  # Save rho, psi and ploidy for future reference
  rho_psi_output = data.frame(rho = c(ascat_optimum_pair$rho,ascat_optimum_pair_fraction_of_genome$rho,ascat_optimum_pair_ref_seg$rho),psi = c(ascat_optimum_pair$psi,ascat_optimum_pair_fraction_of_genome$psi,ascat_optimum_pair_ref_seg$psi), ploidy = c(ascat_optimum_pair$ploidy,ascat_optimum_pair_fraction_of_genome$ploidy,ascat_optimum_pair_ref_seg$ploidy), distance = c(NA,out$distance_without_ref,out$distance), is.best = c(NA,!is.ref.better,is.ref.better),row.names=c("ASCAT","FRAC_GENOME","REF_SEG"))
  write.table(rho_psi_output,paste(outputfile.prefix,"rho_and_psi.txt",sep=""),quote=F,sep="\t")
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
#' @param rho.psi.file String pointing to the rho_and_psi file generated by \code{fit.copy.number}
#' @param output.file Filename of the file where the final copy number fit will be written to
#' @param output.figures.prefix Prefix of the filenames for the chromosome specific copy number figures
#' @param output.gw.figures.prefix Prefix of the filenames for the genome wide copy number figures
#' @param chr_names Vector of allowed chromosome names
#' @param masking_output_file Filename of where the masking details need to be written. Masking is performed to remove very high copy number state segments
#' @param max_allowed_state The maximum CN state allowed (Default 250)
#' @param cn_upper_limit The maximum CN that can be called (Default 1000)
#' @param prior_breakpoints_file A two column file with prior breakpoints, possibly from structural variants. This file must contain two columns: chromosome and position. These are used when making the figures
#' @param gamma Technology specific scaling parameter for LogR (Default 1)
#' @param segmentation.gamma Legacy parameter that is no longer used (Default NA)
#' @param siglevel Threshold under which a p-value becomes significant. When it is significant a second copy number state will be fitted (Default 0.05)
#' @param maxdist Slack in BAF space to allow a segment to be off it's optimum before becoming significant. A segment becomes significant very quickly when a breakpoint is missed, this parameter alleviates the effect (Default 0.01)
#' @param noperms The number of permutations to be run when bootstrapping the confidence intervals on the copy number state of each segment (Default 1000)
#' @param seed Seed to set when performing bootstrapping (Default: Current time)
#' @param calc_seg_baf_option Various options to recalculate the BAF of a segment. Options are: 1 - median, 2 - mean, 3 - ifelse median==0|1, mean, median. (Default: 3)
#' @author dw9, sd11
#' @export

callSubclones = function(sample.name, baf.segmented.file, logr.file, rho.psi.file, output.file, output.figures.prefix, output.gw.figures.prefix, chr_names, masking_output_file, max_allowed_state=250, cn_upper_limit=1000, prior_breakpoints_file=NULL, gamma=1, segmentation.gamma=NA, siglevel=0.05, maxdist=0.01, noperms=1000, seed=as.integer(Sys.time()), calc_seg_baf_option=3) {
  
  set.seed(seed)
  # Load rho/psi/goodness of fit
  res = load.rho.psi.file(rho.psi.file)
  rho = res$rho
  psit = res$psit
  psi = rho*psit + 2 * (1-rho) # psi of all cells
  goodness = res$goodness
  
  # Load the BAF segmented data
  BAFvals = as.data.frame(read_bafsegmented(baf.segmented.file))
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
  LogRvals = as.data.frame(read_logr(logr.file))
  if (colnames(LogRvals)[1] == "X") {
    # If there were rownames, then delete this column. Should not be an issue with new BB runs
    LogRvals = LogRvals[,-1]
  }
  
  # Chromosome names are sometimes 'chr1', etc.
  #if(length(grep("chr",LogRvals[1,1]))>0){
  #	      LogRvals[,1] = gsub("chr","",LogRvals[,1])
  #}
  
  ctrans = c(1:length(chr_names))
  names(ctrans) = chr_names
  ctrans.logR = c(1:length(chr_names))
  names(ctrans.logR) = chr_names
  
  # = as.vector(ctrans.logR[as.vector(LogRvals[,1])]*1000000000+LogRvals[,2])
  BAFpos = as.vector(ctrans[as.vector(BAFvals[,1])]*1000000000+BAFvals[,2])
  
  ################################################################################################
  # Determine copy number for each segment
  ################################################################################################
  res = determine_copynumber(BAFvals, LogRvals, rho, psi, gamma, ctrans, ctrans.logR, maxdist, siglevel, noperms, cn_upper_limit)
  subcloneres = res$subcloneres
  #write.table(subcloneres, gsub(".txt", "_1.txt", output.file), quote=F, col.names=T, row.names=F, sep="\t")
  write.table(subcloneres, paste0(tools::file_path_sans_ext(output.file),"_1.",tools::file_ext(output.file),sep=""), quote=F, col.names=T, row.names=F, sep="\t")
  # Scan the segments for cases that should be merged
  res = merge_segments(subcloneres, BAFvals, LogRvals, rho, psi, gamma, calc_seg_baf_option)
  BAFvals = res$bafsegmented
  
  res = determine_copynumber(BAFvals, LogRvals, rho, psi, gamma, ctrans, ctrans.logR, maxdist, siglevel, noperms, cn_upper_limit)
  subcloneres = res$subcloneres
  BAFpvals = res$BAFpvals
  
  # Scan for very high copy number segments and set those to NA - This is in part an artifact of small segments
  res = mask_high_cn_segments(subcloneres, BAFvals, max_allowed_state)
  subcloneres = res$subclones
  # No longer writing out the BAFsegmented data after masking
  #BAFvals = res$bafsegmented
  #write.table(BAFvals, file=baf.segmented.file, sep="\t", row.names=F, col.names=T, quote=F)
  # Write the masking details to file
  masking_details = data.frame(samplename=sample.name, masked_count=res$masked_count, masked_size=res$masked_size, max_allowed_state=max_allowed_state)
  write.table(masking_details, file=masking_output_file, quote=F, col.names=T, row.names=F, sep="\t")
  
  # Write the final copy number profile 
  # NAP: generating two output files: first reporting solution A and the second reporting alternative solutions (B to F)
  write.table(subcloneres[,c(1:3,8:13)], output.file, quote=F, col.names=T, row.names=F, sep="\t")

  #write.table(subcloneres, gsub(".txt","_extended.txt",output.file), quote=F, col.names=T, row.names=F, sep="\t")
  write.table(subcloneres, paste0(tools::file_path_sans_ext(output.file),"_extended.",tools::file_ext(output.file),sep=""), quote=F, col.names=T, row.names=F, sep="\t")

  # NAP - November 2023
  # Recalculate PGA.is.clonal to match the final copy number profile in copynumber.txt file (previously subclones.txt file)
  subcloneres$length = subcloneres$endpos-subcloneres$startpos
  subcloneres_subclonal = subcloneres[which(subcloneres$frac1_A<1),]
  diploid = which(subcloneres$nMaj1_A==1 & subcloneres$nMin1_A==1 & subcloneres$frac1_A==1)
  # NAP - June 2025 
  # Check 'diploid' length for rare edge cases
  if (length(diploid) > 0) {
    cna = subcloneres[-diploid,]
  } else {
    cna = subcloneres
    print("No diploid region found in copy number profile - likely due to WGD or error in fitting copy number in rare cases") 
  }
  
  if(nrow(cna) == 0 || sum(cna$length) == 0) {
    # No copy number alterations found
    goodness <- 1.0  # 100% clonal (no CNAs to be subclonal)
    print("No copy number alterations detected - setting PGA.is.clonal to 100%\n")
  } else if(nrow(subcloneres_subclonal) == 0) {
    # No subclonal segments
    goodness <- 1.0  # 100% clonal
    print("No subclonal segments detected - setting PGA.is.clonal to 100%\n")
  } else {
    subclonal_fraction <- sum(subcloneres_subclonal$length) / sum(cna$length)
    goodness <- 1 - subclonal_fraction
  
    # Ensure goodness is within valid range [0,1]
    goodness <- max(0, min(1, goodness))
  }
  print(paste0("PGA.is.clonal = ",sprintf("%2.1f",goodness*100),"%"))
  
  ################################################################################################
  # Make a plot per chromosome
  ################################################################################################
  # Collapse the BAFsegmented into breakpoints to be used in plotting
  segment_breakpoints = collapse_bafsegmented_to_segments(BAFvals)
  if (!is.null(prior_breakpoints_file) & !ifelse(is.null(prior_breakpoints_file), TRUE, prior_breakpoints_file=="NA") & !ifelse(is.null(prior_breakpoints_file), TRUE, is.na(prior_breakpoints_file))) {
    svs = read.table(prior_breakpoints_file, header=T, stringsAsFactors=F)
  }
  
  # Create a plot per chromosome that shows the segments with their CN state in text
  for (chr in chr_names) {
    pos = SNPpos[SNPpos[,1]==chr, 2]
    #if no points to plot, skip
    if (length(pos)==0) { next }
    
    if (!is.null(prior_breakpoints_file) & !ifelse(is.null(prior_breakpoints_file), TRUE, prior_breakpoints_file=="NA") & !ifelse(is.null(prior_breakpoints_file), TRUE, is.na(prior_breakpoints_file))) {
      svs_pos = svs[svs$chromosome==chr,]$position / 1000000
    } else {
      svs_pos = NULL
    }
    
    breakpoints_pos = segment_breakpoints[segment_breakpoints$chromosome==chr,]
    breakpoints_pos = sort(unique(c(breakpoints_pos$start, breakpoints_pos$end) / 1000000))
    
    png(filename = paste(output.figures.prefix, chr,".png",sep=""), width = 2000, height = 2000, res = 200, type = "cairo")
    create.subclonal.cn.plot(chrom=chr,
                             chrom.position=pos/1000000,
                             LogRposke=LogRvals[LogRvals[,1]==chr,2],
                             LogRchr=LogRvals[LogRvals[,1]==chr,3],
                             BAFchr=BAF[SNPpos[,1]==chr],
                             BAFsegchr=BAFseg[SNPpos[,1]==chr],
                             BAFpvalschr=BAFpvals[SNPpos[,1]==chr],
                             subcloneres=subcloneres,
                             breakpoints_pos=breakpoints_pos,
                             svs_pos=svs_pos,
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
  
  # Recalculate the ploidy based on the actual fit
  seg_length = floor((subclones$endpos-subclones$startpos)/1000)
  is_subclonal_maj = abs(subclones$nMaj1_A - subclones$nMaj2_A) > 0
  is_subclonal_min = abs(subclones$nMin1_A - subclones$nMin2_A) > 0
  is_subclonal_maj[is.na(is_subclonal_maj)] = F
  is_subclonal_min[is.na(is_subclonal_min)] = F
  segment_states_min = subclones$nMin1_A * ifelse(is_subclonal_min, subclones$frac1_A, 1)  + ifelse(is_subclonal_min, subclones$nMin2_A, 0) * ifelse(is_subclonal_min, subclones$frac2_A, 0)
  segment_states_maj = subclones$nMaj1_A * ifelse(is_subclonal_maj, subclones$frac1_A, 1)  + ifelse(is_subclonal_maj, subclones$nMaj2_A, 0) * ifelse(is_subclonal_maj, subclones$frac2_A, 0)
  ploidy = sum((segment_states_min+segment_states_maj) * seg_length, na.rm=T) / sum(seg_length, na.rm=T)
  
  # Plot genome wide figures
  plot.gw.subclonal.cn(subclones=subclones, BAFvals=BAFvals, rho=rho, ploidy=ploidy, goodness=goodness, output.gw.figures.prefix=output.gw.figures.prefix, chr.names=chr_names, tumourname=sample.name)
  
  # Create user friendly cellularity and ploidy output file
  cellularity_ploidy_output = data.frame(purity = c(rho), ploidy = c(ploidy), psi = c(psit))

  # cellularity_file = gsub("_.+\\.txt$", "_purity_ploidy.txt", output.file) # NAP: updated the name of the output file, consistent with new title (and added flexibility with what output.file is named)
  cellularity_file = paste0(sample.name,"_purity_ploidy.txt") 

  write.table(cellularity_ploidy_output, cellularity_file, quote=F, sep="\t", row.names=F)
}


#' Given all the determined values make a copy number call for each segment
#'
#' @param BAFvals BAFsegmented data.frame with 5 columns
#' @param LogRvals Raw logR values in data.frame with 3 columns
#' @param rho Optimal rho value, the choosen cellularity
#' @param psi Optimal psi value, the choosen ploidy
#' @param gamma Platform gamma parameter
#' @param ctrans Named vector of chromosome names
#' @param ctrans.logR Named vector of chromosome names
#' @param maxdist Max distance a segment is tolerated to be not considered for subclonal copy number
#' @param siglevel Level at which a segment can become significantly different from the nearest clonal state
#' @param noperms Number of bootstrap permutations
#' @param cn_upper_limit Maximum number of CN that can be called
#' @return A data.frame with copy number determined for each segment
#' @author dw9
#' @noRd
determine_copynumber = function(BAFvals, LogRvals, rho, psi, gamma, ctrans, ctrans.logR, maxdist, siglevel, noperms, cn_upper_limit) {
  BAFphased = BAFvals[,4]
  BAFseg = BAFvals[,5]
  BAFpos = as.vector(ctrans[as.vector(BAFvals[,1])]*1000000000+BAFvals[,2])
  LogRpos = as.vector(ctrans.logR[as.vector(LogRvals[,1])]*1000000000+LogRvals[,2])
  
  #DCW 240314
  switchpoints = c(0,which(BAFseg[-1] != BAFseg[-(length(BAFseg))] | BAFvals[-1,1] != BAFvals[-nrow(BAFvals),1]),length(BAFseg))
  BAFlevels = BAFseg[switchpoints[-1]]
  
  pval = NULL
  BAFpvals = vector(length=length(BAFseg))
  subcloneres = NULL
  
  for (i in 1:length(BAFlevels)) {
    # subcloneres = rbind(subcloneres, fit_segment(BAFpos, LogRpos, BAFlevels, BAFphased, LogRvals, switchpoints, rho, psi, gamma, i))
    l = BAFlevels[i]
    
    # Make sure that BAF>=0.5, otherwise nMajor and nMinor may be the wrong way around
    l = max(l,1-l)
    
    BAFke = BAFphased[(switchpoints[i]+1):switchpoints[i+1]]
    
    #startpos = min(BAFpos[names(BAFke)])
    #endpos = max(BAFpos[names(BAFke)])
    startpos = min(BAFpos[(switchpoints[i]+1):switchpoints[i+1]])
    endpos = max(BAFpos[(switchpoints[i]+1):switchpoints[i+1]])
    #chrom = names(ctrans[floor(startpos/1000000000)])
    # Assuming all SNPs in this segment are on the same chromosome
    chrom = BAFvals[(switchpoints[i]+1):switchpoints[i+1],]$Chromosome[1]
    LogR = mean(LogRvals[LogRpos>=startpos&LogRpos<=endpos & !is.infinite(LogRvals[,3]),3],na.rm=T)
    
    # if we don't have a value for LogR, fill in 0
    if (is.na(LogR)) {
      LogR = 0
    }
    nMajor = (rho-1+l*psi*2^(LogR/gamma))/rho
    nMinor = (rho-1+(1-l)*psi*2^(LogR/gamma))/rho
   
    # Occasionally nMinor can be NA due to zero coverage, skip when this occurs
    if (is.na(nMinor)) {
      next
    }

    # Increase nMajor and nMinor together, to avoid impossible combinations (with negative subclonal fractions)
    if (nMinor<0) {
      if (l==1) {
        # Avoid calling infinite copy number
        nMajor = cn_upper_limit
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
    if (is.na(sd(BAFke)) || sd(BAFke)==0) {
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
      
      subcloneres = rbind(subcloneres, c(chrom,startpos-floor(startpos/1000000000)*1000000000,
                                         endpos-floor(endpos/1000000000)*1000000000,l,pval[i],LogR,ntot,
                                         nMaj1[1],nMin1[1],tau[1],nMaj2[1],nMin2[1],1-tau[1],sdtau[1],sdtaubootstrap[1],tau25[1],tau975[1],
                                         nMaj1[2],nMin1[2],tau[2],nMaj2[2],nMin2[2],1-tau[2],sdtau[2],sdtaubootstrap[2],tau25[2],tau975[2],
                                         nMaj1[3],nMin1[3],tau[3],nMaj2[3],nMin2[3],1-tau[3],sdtau[3],sdtaubootstrap[3],tau25[3],tau975[3],
                                         nMaj1[4],nMin1[4],tau[4],nMaj2[4],nMin2[4],1-tau[4],sdtau[4],sdtaubootstrap[4],tau25[4],tau975[4],
                                         nMaj1[5],nMin1[5],tau[5],nMaj2[5],nMin2[5],1-tau[5],sdtau[5],sdtaubootstrap[5],tau25[5],tau975[5],
                                         nMaj1[6],nMin1[6],tau[6],nMaj2[6],nMin2[6],1-tau[6],sdtau[6],sdtaubootstrap[6],tau25[6],tau975[6]))
    }else {
      #if called as clonal, use the best corner from the nearest edge
      subcloneres = rbind(subcloneres, c(chrom,startpos-floor(startpos/1000000000)*1000000000,
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
  subcloneres = as.data.frame(subcloneres)
  for (i in 2:ncol(subcloneres)) {
    subcloneres[,i] = as.numeric(as.character(subcloneres[,i]))
  }
  return(list(subcloneres=subcloneres, BAFpvals=BAFpvals))
}


#' Merge copy number segments
#'
#' Merges segments if there is not enough evidence for them to be separate. Two adjacent segments are merged
#' when they are either fit with the same clonal copy number state or when their BAF is not significantly different
#' and their logR puts them in the same square.
#' @param subclones A completely fit copy number profile in Battenberg output format
#' @param bafsegmented A BAFsegmented data.frame with the 5 columns that corresponds to the subclones file
#' @param logR The raw logR data
#' @param rho The rho estimate that the profile was fit with
#' @param psi the psi estimate that the profile was fit with
#' @param platform_gamma The gamma parameter for this platform
#' @param calc_seg_baf_option Various options to recalculate the BAF of a segment. Options are: 1 - median, 2 - mean, 3 - ifelse median== 0|1, mean, median. (Default: 3)
#' @param verbose A boolean to show merging operations (Default: FALSE)
#' @return A list with two fields: bafsegmented and subclones. The subclones field contains a data.frame in
#' Battenberg output format with the merged segments. The bafsegmented field contains the BAFsegmented data
#' corresponding to the provided subclones data.frame.
#' @author sd11, tl
#' @noRd
merge_segments=function(subclones, bafsegmented, logR, rho, psi, platform_gamma, calc_seg_baf_option=3, verbose=F) {
  calc_nmin = function(rho, psi, baf, logr, platform_gamma) {
    return((rho-1-(baf-1)*2^(logr/platform_gamma)*((1-rho)*2+rho*psi))/rho)
  }
  calc_nmaj = function(rho, psi, baf, logr, platform_gamma) {
    return((rho-1+baf*2^(logr/platform_gamma)*((1-rho)*2+rho*psi))/rho)
  }
  # Convert DF into GRanges objects
  df2gr=function(DF,chr,pos1,pos2) {
    return(GenomicRanges::makeGRangesFromDataFrame(df=DF,
                                                   keep.extra.columns=T,
                                                   ignore.strand=T,
                                                   seqinfo=NULL,
                                                   seqnames.field=chr,
                                                   start.field=pos1,
                                                   end.field=pos2,
                                                   starts.in.df.are.0based=F))
  }
  # Function called when two segments have not been merged so there is no need to recheck those again
  updateNeighbour=function(subclones,INDEX,INDEX_N) {
    if (INDEX_N>INDEX) {
      subclones$Next_checked[INDEX]=T
      subclones$Prev_checked[INDEX_N]=T
    } else {
      subclones$Prev_checked[INDEX]=T
      subclones$Next_checked[INDEX_N]=T
    }
    return(subclones)
  }
  # Function called when two segments have been merged so we need to recheck its two neighbours
  updateAround=function(subclones,INDEX) {
    if (INDEX>1) {
      subclones$Prev_checked[INDEX]=F
      subclones$Next_checked[INDEX-1]=F
    } else {
      subclones$Prev_checked[INDEX]=T
    }
    if (INDEX<length(subclones)) {
      subclones$Next_checked[INDEX]=F
      subclones$Prev_checked[INDEX+1]=F
    } else {
      subclones$Next_checked[INDEX]=T
    }
    return(subclones)
  }
  # Function called to test whether two segments must be checked
  checkStatus=function(subclones,INDEX,INDEX_N) {
    if (INDEX_N>INDEX) {
      # Largest segment (INDEX_N) is after smallest one (INDEX)
      stopifnot(subclones$Next_checked[INDEX]==subclones$Prev_checked[INDEX_N])
      if (subclones$Next_checked[INDEX] && subclones$Prev_checked[INDEX_N]) {
        return(T)
      } else {
        return(F)
      }
    } else {
      # Largest segment (INDEX_N) is before smallest one (INDEX)
      stopifnot(subclones$Prev_checked[INDEX]==subclones$Next_checked[INDEX_N])
      if (subclones$Prev_checked[INDEX] && subclones$Next_checked[INDEX_N]) {
        return(T)
      } else {
        return(F)
      }
    }
  }
  # Function to merge two segments
  merge_seg=function(subclones,bafsegmented,logR,INDEX,INDEX_N,calc_seg_baf_option) {
    # Update start/end information
    if (INDEX_N<INDEX) {
      GenomicRanges::end(subclones[INDEX_N])=GenomicRanges::end(subclones[INDEX])
    } else {
      GenomicRanges::start(subclones[INDEX_N])=GenomicRanges::start(subclones[INDEX])
    }
    # Remove segment
    subclones=subclones[-INDEX]
    if(INDEX_N<INDEX) INDEX=INDEX-1
    # Reset neighbour checking
    subclones=updateAround(subclones,INDEX)
    if (calc_seg_baf_option==1) {
      # This uses median
      NEW_BAF = median(bafsegmented$BAFphased[GenomicRanges::findOverlaps(subclones[INDEX],bafsegmented)@to], na.rm=T)
    } else if (calc_seg_baf_option==2) {
      # This uses mean
      NEW_BAF = mean(bafsegmented$BAFphased[GenomicRanges::findOverlaps(subclones[INDEX],bafsegmented)@to], na.rm=T)
    } else if (calc_seg_baf_option==3) {
      # We'll prefer the median BAF as a segment summary
      # but change to the mean when the median is extreme
      # as at 0 or 1 the BAF is uninformative for the fitting
      median_BAF = median(bafsegmented$BAFphased[GenomicRanges::findOverlaps(subclones[INDEX],bafsegmented)@to], na.rm=T)
      mean_BAF = mean(bafsegmented$BAFphased[GenomicRanges::findOverlaps(subclones[INDEX],bafsegmented)@to], na.rm=T)
      if (median_BAF!=0 && median_BAF!=1) {
        NEW_BAF = median_BAF
      } else {
        NEW_BAF = mean_BAF
      }
      rm(median_BAF,mean_BAF)
    }
    # Update both BAF and logR information
    subclones[INDEX]$BAF=NEW_BAF
    INDEX_logR=GenomicRanges::findOverlaps(subclones[INDEX],logR)@to
    if (length(INDEX_logR)==0) {
      subclones[INDEX]$LogR=0
    } else {
      subclones[INDEX]$LogR=mean(logR$logR[INDEX_logR], na.rm=T)
    }
    rm(INDEX_logR)
    # Update segmented baf
    bafsegmented$BAFseg[GenomicRanges::findOverlaps(subclones[INDEX],bafsegmented)@to] = NEW_BAF
    # TODO Update the logRseg as well
    # Reset IDs
    subclones$ID=1:length(subclones)
    return(list(subclones=subclones,bafsegmented=bafsegmented))
  }
  requireNamespace("GenomicRanges")
  if (!(calc_seg_baf_option %in% 1:3)) calc_seg_baf_option=3
  # Convert DFs into GRanges objects
  if (verbose) print('Convert DFs into GRanges objects')
  subclones=df2gr(subclones,'chr','startpos','endpos')
  str(bafsegmented)
  bafsegmented=df2gr(bafsegmented,'Chromosome','Position','Position')
  str(logR)
  logR=df2gr(logR,'Chromosome','Position','Position')
  names(GenomicRanges::mcols(logR))='logR'
  # Split GRanges objects by chromosomes
  chr_names=GenomicRanges::seqnames(GenomicRanges::seqinfo(bafsegmented))
  subclones=lapply(chr_names,function(x) subclones[GenomicRanges::seqnames(subclones)==x])
  bafsegmented=lapply(chr_names,function(x) bafsegmented[GenomicRanges::seqnames(bafsegmented)==x])
  logR=lapply(chr_names,function(x) logR[GenomicRanges::seqnames(logR)==x])
  
  stopifnot(all(sapply(subclones,length)>0) && all(sapply(bafsegmented,length)>0) && all(sapply(logR,length)>0))
  names(subclones)=chr_names
  names(bafsegmented)=chr_names
  names(logR)=chr_names
  # For each chromosome
  for (CHR in chr_names) {
    if (verbose) print(paste0('Merging segments within: ',CHR))
    # Define ID, Prev_checked and Next_checked to help processing data 
    subclones[[CHR]]$ID=1:length(subclones[[CHR]])
    subclones[[CHR]]$Prev_checked=F
    subclones[[CHR]]$Next_checked=F
    subclones[[CHR]]$Prev_checked[1]=T
    subclones[[CHR]]$Next_checked[length(subclones[[CHR]])]=T
    # Pick all possible IDs
    IDs=subclones[[CHR]]$ID
    while (length(IDs)!=0) {
      # Amongst all IDs, select the ones that must be checked
      IDs=subclones[[CHR]]$ID[which(!subclones[[CHR]]$Prev_checked | !subclones[[CHR]]$Next_checked)]
      if (length(IDs)==0) break
      # Amongst all of those, select the smallest one
      INDEX=IDs[which.min(GenomicRanges::width(subclones[[CHR]][which(subclones[[CHR]]$ID %in% IDs)]))]
      # Select neighbours (two or one if segments is first or last)
      if (INDEX==1) {
        Neighbours=order(GenomicRanges::distance(subclones[[CHR]][INDEX],subclones[[CHR]][INDEX+1]))
        names(Neighbours)=INDEX+1
      } else if (INDEX==length(subclones[[CHR]])) {
        Neighbours=order(GenomicRanges::distance(subclones[[CHR]][INDEX],subclones[[CHR]][INDEX-1]))
        names(Neighbours)=INDEX-1
      } else {
        Neighbours=order(GenomicRanges::distance(subclones[[CHR]][INDEX],subclones[[CHR]][INDEX+c(-1,1)]))
        names(Neighbours)=INDEX+c(-1,1)
      }
      if (verbose) print(paste0('Working on segment: ',INDEX,' (',subclones[[CHR]][INDEX],')'))
      # For each neighbour
      for (i in Neighbours) {
        INDEX_N=as.numeric(names(Neighbours[i]))
        if (verbose) print(paste0('Checking neighbour: ',INDEX_N,' (',subclones[[CHR]][INDEX_N],'; distance=',GenomicRanges::distance(subclones[[CHR]][INDEX],subclones[[CHR]][INDEX_N]),')'))
        # Test whether seg and neighbour (INDEX and INDEX_N) have already been checked
        if (checkStatus(subclones[[CHR]],INDEX,INDEX_N)) {if (verbose) {print('Already checked')}; next}
        # Test whether seg and neighbour are far away from each other
        if (GenomicRanges::distance(subclones[[CHR]][INDEX],subclones[[CHR]][INDEX_N])>3e6) {
          if (verbose) print('Distance > 3Mb - do not merge')
          subclones[[CHR]]=updateNeighbour(subclones[[CHR]],INDEX,INDEX_N)
        } else {
          # Test whether seg and neighbour have the same clonal CN solution
          if (subclones[[CHR]]$nMaj1_A[INDEX]==subclones[[CHR]]$nMaj1_A[INDEX_N] && subclones[[CHR]]$nMin1_A[INDEX]==subclones[[CHR]]$nMin1_A[INDEX_N] && subclones[[CHR]]$frac1_A[INDEX]==1 && subclones[[CHR]]$frac1_A[INDEX_N]==1) {
            if (verbose) print('Same clonal CN solution - merge')
            res=merge_seg(subclones[[CHR]],bafsegmented[[CHR]],logR[[CHR]],INDEX,INDEX_N,calc_seg_baf_option)
            subclones[[CHR]]=res$subclones
            bafsegmented[[CHR]]=res$bafsegmented
            rm(res)
            break
          } else {
            # Test whether seg and neighbour have different BAF/logR distributions
            if (verbose) print('Different CN solutions: check BAF and logR')
            nmin_curr = round(calc_nmin(rho, psi, subclones[[CHR]]$BAF[INDEX], subclones[[CHR]]$LogR[INDEX], platform_gamma))
            nmaj_curr = round(calc_nmaj(rho, psi, subclones[[CHR]]$BAF[INDEX], subclones[[CHR]]$LogR[INDEX], platform_gamma))
            nmin_other = round(calc_nmin(rho, psi, subclones[[CHR]]$BAF[INDEX_N], subclones[[CHR]]$LogR[INDEX_N], platform_gamma))
            nmaj_other = round(calc_nmaj(rho, psi, subclones[[CHR]]$BAF[INDEX_N], subclones[[CHR]]$LogR[INDEX_N], platform_gamma))
            if (nmin_curr==nmin_other || nmaj_curr==nmaj_other) {
              # Test whether there are more than 10 values to check significance
              if (sum(!is.na(logR[[CHR]]$logR[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX],logR[[CHR]])@to])) > 10 &&
                  sum(!is.na(logR[[CHR]]$logR[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX_N],logR[[CHR]])@to])) > 10 &&
                  sum(!is.na(bafsegmented[[CHR]]$BAFphased[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX],bafsegmented[[CHR]])@to])) > 10 &&
                  sum(!is.na(bafsegmented[[CHR]]$BAFphased[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX_N],bafsegmented[[CHR]])@to])) > 10) {
                logr_significant = t.test(logR[[CHR]]$logR[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX],logR[[CHR]])@to],
                                          logR[[CHR]]$logR[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX_N],logR[[CHR]])@to])$p.value < 0.05
                baf_significant = t.test(bafsegmented[[CHR]]$BAFphased[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX],bafsegmented[[CHR]])@to],
                                         bafsegmented[[CHR]]$BAFphased[GenomicRanges::findOverlaps(subclones[[CHR]][INDEX_N],bafsegmented[[CHR]])@to])$p.value < 0.05
                if ((!logr_significant) && (!baf_significant)) {
                  if (verbose) print('No significant difference - merge')
                  res=merge_seg(subclones[[CHR]],bafsegmented[[CHR]],logR[[CHR]],INDEX,INDEX_N,calc_seg_baf_option)
                  subclones[[CHR]]=res$subclones
                  bafsegmented[[CHR]]=res$bafsegmented
                  rm(res)
                  break
                } else {
                  if (verbose) print('Significant difference - do not merge')
                  subclones[[CHR]]=updateNeighbour(subclones[[CHR]],INDEX,INDEX_N)
                }
              } else {
                if (verbose) print('Too few values - do not merge')
                subclones[[CHR]]=updateNeighbour(subclones[[CHR]],INDEX,INDEX_N)
              }
            } else {
              if (verbose) print('Different squares - do not merge')
              subclones[[CHR]]=updateNeighbour(subclones[[CHR]],INDEX,INDEX_N)
            }
          }
        }
      }; rm(i)
    }
  }; rm(CHR)
  if (verbose) print('Convert GRanges objects into DFs')
  bafsegmented=data.frame(Reduce(c,bafsegmented),stringsAsFactors=F)[,-c(3:5)]
  bafsegmented$seqnames=as.character(bafsegmented$seqnames)
  colnames(bafsegmented)[1:2]=c('Chromosome','Position')
  subclones=data.frame(Reduce(c,subclones),stringsAsFactors=F)[,-c(4:5)]
  subclones$seqnames=as.character(subclones$seqnames)
  colnames(subclones)[1:3]=c('chr','startpos','endpos')
  subclones$ID=NULL
  subclones$Prev_checked=NULL
  subclones$Next_checked=NULL
  return(list(bafsegmented=bafsegmented, subclones=subclones))
}

#' Mask segments that have a too high CN state
#' @param subclones Subclones output data
#' @param bafsegmented BAFsegmented data
#' @param max_allowed_state The maximum state allowed before overruling takes place
#' @return A list with the masked subclones, bafsegmented and the number of segments masked and their total genome size
#' @author sd11
mask_high_cn_segments = function(subclones, bafsegmented, max_allowed_state) {
  count = 0
  masked_size = 0
  for (i in 1:nrow(subclones)) {
    if (subclones$nMaj1_A[i] > max_allowed_state | subclones$nMin1_A[i] > max_allowed_state) {
      # Mask this segment
      subclones[i, "nMaj1_A"] = NA
      subclones[i, "nMin1_A"] = NA
      subclones[i, "nMaj2_A"] = NA
      subclones[i, "nMin2_A"] = NA
      # Mask the BAFsegmented
      bafsegmented[subclones$chr[i] == bafsegmented$Chromosome & subclones$startpos[i] < bafsegmented$Position & subclones$endpos[i] >= bafsegmented$Position,c("BAFseg")] = NA
      count = count+1
      masked_size = masked_size + (subclones$endpos[i]-subclones$startpos[i])
    }
  }
  return(list(subclones=subclones, bafsegmented=bafsegmented, masked_count=count, masked_size=masked_size))
}


#' Plot the copy number genome wide in two different ways. This creates the Battenberg average
#' profile where subclonal copy number is represented as a mixture of two different states and
#' the Battenberg subclones profile where subclonal copy number is plotted as two different
#' separate states. The thickness of the line represents the fraction of tumour cells carying
#' the particular state.
#' @noRd
plot.gw.subclonal.cn = function(subclones, BAFvals, rho, ploidy, goodness, output.gw.figures.prefix, chr.names, tumourname) {
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
  png(filename = paste(output.gw.figures.prefix, "_average.png", sep=""), width = 2000, height = 500, res = 200, type = "cairo")
  create.bb.plot.average(bafsegmented=BAFvals,
                         ploidy=ploidy,
                         rho=rho,
                         goodnessOfFit=goodness,
                         pos_min=pos_min,
                         pos_max=pos_max,
                         segment_states_min=segment_states_min,
                         segment_states_tot=segment_states_tot,
                         chr.segs=chr.segs,
                         chr.names=chr.names,
                         tumourname=tumourname)
  dev.off()
  
  # Plot subclonal copy number as two separate states
  png(filename = paste(output.gw.figures.prefix, "_subclones.png", sep=""), width = 2000, height = 500, res = 200, type = "cairo")
  create.bb.plot.subclones(bafsegmented=BAFvals,
                           subclones=subclones,
                           ploidy=ploidy,
                           rho=rho,
                           goodnessOfFit=goodness,
                           pos_min=pos_min,
                           pos_max=pos_max,
                           subcl_min=subcl_min,
                           subcl_max=subcl_max,
                           is_subclonal=is_subclonal,
                           is_subclonal_maj=is_subclonal_maj,
                           is_subclonal_min=is_subclonal_min,
                           chr.segs=chr.segs,
                           chr.names=chr.names,
                           tumourname=tumourname)
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

#' Collapse a BAFsegmented file into segment start and end points
#'
#' This function looks through the BAFsegmented for stretches of equal
#' BAFseg and records the start and end coordinates in a data.frame
#' @param bafsegmented The BAFsegmented output from segmentation
#' @return A data.frame with columns chromosome, start and end
#' @author sd11
#' @noRd
collapse_bafsegmented_to_segments = function(bafsegmented) {
  segments_collapsed = data.frame()
  for (chrom in unique(bafsegmented$Chromosome)) {
    bafsegmented_chrom = bafsegmented[bafsegmented$Chromosome==chrom,]
    segments = rle(bafsegmented_chrom$BAFseg)
    startpoint = 1
    for (i in 1:length(segments$lengths)) {
      endpoint = startpoint+segments$lengths[i]-1
      segments_collapsed = rbind(segments_collapsed,
                                 data.frame(chromosome=chrom, start=bafsegmented_chrom$Position[startpoint], end=bafsegmented_chrom$Position[endpoint]))
      startpoint = endpoint+1
    }
  }
  return(segments_collapsed)
}

#' Function to make additional figures
#'
#' @param samplename Name of the sample for the plot title
#' @param logr_file File containing all logR data
#' @param bafsegmented_file File containing the BAFsegmented data
#' @param logrsegmented_file File with the logRsegmented data
#' @param allelecounts_file Optional file with raw allele counts (Default: NULL)
#' @author sd11
#' @export
make_posthoc_plots = function(samplename, logr_file, bafsegmented_file, logrsegmented_file, allelecounts_file=NULL) {
  # Make some post-hoc plots
  logr = Battenberg::read_table_generic(logr_file)
  bafsegmented = as.data.frame(Battenberg::read_table_generic(bafsegmented_file))
  logrsegmented = as.data.frame(Battenberg::read_table_generic(logrsegmented_file, header=F))
  colnames(logrsegmented) = c("Chromosome", "Position", "logRseg")
  outputfile = paste0(samplename, "_alleleratio.png")
  allele_ratio_plot(samplename=samplename, logr=logr, bafsegmented=bafsegmented, logrsegmented=logrsegmented, outputfile=outputfile, max.plot.cn=8)
  
  if (!is.null(allelecounts_file)) {
    allelecounts = as.data.frame(Battenberg::read_table_generic(allelecounts_file))
    outputfile = paste0(samplename, "_coverage.png")
    coverage_plot(samplename, allelecounts, outputfile)
  }
}


#' Fit ChrX subclonal copy number (male only)
#'
#' Function to call ChrX copy number based on LogR (suitable for male samples). Copy number
#' cannot be called for the non-PAR region of ChrX due to the hemizygosity of all 1000G SNPs.
#' This function enables calling subclonal copy number for the non-PAR region by segmenting LogR.
#' A number of correction steps are undertaken to account for the noisy nature of LogR. This function
#' requires the following libraries: copynumber, data.table and ggplot2. It reads in three files generated
#' by previous steps of Battenberg, namely samplename_mutantLogR_gcCorrected.tab, samplename_purity_ploidy.txt
#' and samplename_copynumber_extended.txt.
#' This function will also update the Battenberg genome-wide profile plots (average.png and subclones.png) to include the chrX profile by also
#' reading in the samplename.BAFsegmented.txt and samplename_rho_psi.txt files
#' @param tumourname The sample name used for Battenberg (i.e. the tumour BAM file name without the .bam extension)
#' @param X_gamma The PCF gamma value for segmentation of 1000G SNP LogR values (Default 1000)
#' @param X_kmin The min number of SNPs to support a segment in PCF of LogR values (Default 100)
#' @param genomebuild The genome build used in running Battenberg (hg19 or hg38)
#' @param AR Should the segment carrying the androgen receptor (AR) locus to be visually distinguished in average plot? (Default TRUE)
#' @param prior_breakpoints_file A two column text file with prior genome-wide breakpoints, possibly from structural variants. This file must contain two columns with headers "chr" and "pos" representing chromosome and position.
#' @param chrom_names A vector containing the names of chromosomes to be included in the final genome-wide Battenberg copy number plot with chrX
#' @author naser.ansari-pour
#' @export

callChrXsubclones = function(tumourname,X_gamma=1000,X_kmin=100,genomebuild,AR=TRUE,prior_breakpoints_file=NULL,chrom_names,data_type="wgs"){
  
  print(tumourname)
  
  if (genomebuild=="hg19"){
    par_regions=c(2699520,155260560)
    x_centromere=c(58632012,61632012)
    ar=data.frame(startpos=66763874,endpos=66950461)
  } else if (genomebuild=="hg38") {
    par_regions=c(2781479,156030895)
    x_centromere=c(58605580,62412542)
    ar=data.frame(startpos=67544021,endpos=67730619)
  } else {
    stop("Genomebuild not supported for callChrXsubclones")
  }
  
  if (data_type=="wgs" | data_type=="WGS") {
    PCFinput=data.frame(read_table_generic(paste0(tumourname,"_mutantLogR_gcCorrected.tab")),stringsAsFactors=F)
  } else {
    PCFinput=data.frame(read_table_generic(paste0(tumourname,"_mutantLogR.tab")),stringsAsFactors=F)
  }
  ChrNotation=unique(PCFinput[which(!is.na(match(PCFinput$Chromosome,c("X","chrX")))),]$Chromosome) # find the chromosome notation
  PCFinput=PCFinput[which(PCFinput$Chromosome==ChrNotation & PCFinput$Position>par_regions[1] & PCFinput$Position<par_regions[2]),] # get nonPAR using par_regions based on genomebuild
  colnames(PCFinput)[3]=tumourname
  print(paste("Number of chrX nonPAR SNPs =",nrow(PCFinput)))
  
  if (!is.null(prior_breakpoints_file)) {
    sv=read.table(prior_breakpoints_file, header=T, stringsAsFactors=F)
    sv=sv[which(!is.na(match(sv$chr,c("X","chrX")))),]
    # check if there are breakpoints within chrX
    if (nrow(sv)>0){
      # make sure all SV breakpoint positions are within the LogR data range and not outside of it  
      svpos=sv[which((sv$pos > min(PCFinput$Position)) & (sv$pos < max(PCFinput$Position))),"pos"]
      breakpoints=c(min(PCFinput$Position),svpos,max(PCFinput$Position))
      PCF=data.frame()
      for (j in 1:(length(breakpoints)-1)) {
        PCFinput_sv=PCFinput[which(PCFinput$Position>=breakpoints[j] & PCFinput$Position<breakpoints[j+1]),]
        # in case there is no SNP between two SVs on chrX
        if (nrow(PCFinput_sv)==0) next
        PCF_sv=copynumber::pcf(PCFinput_sv,gamma=X_gamma,kmin=X_kmin)
        PCF=rbind(PCF,PCF_sv)
      }
    } else {
      PCF=copynumber::pcf(PCFinput,gamma=X_gamma,kmin=X_kmin)
    }
  } else {
    PCF=copynumber::pcf(PCFinput,gamma=X_gamma,kmin=X_kmin)
  }
  write.table(PCF,paste0(tumourname,"_PCF_gamma_",X_gamma,"_chrX.txt"),col.names=T,row.names=F,quote=F,sep="\t")
  print("PCF segmentation done")
  
  # INPUT for copy number inference
  SAMPLEsegs=data.frame(PCF,stringsAsFactors=F)
  pupl=read.table(paste0(tumourname,"_purity_ploidy.txt"),header=T,stringsAsFactors=F)
  SAMPLEpurity=pupl[,1] # SAMPLEpurity=pupl$cellularity in previous Battenberg version; change from pupl$purity to pupl[,1] for universality
  #SAMPLEwgd=ifelse(round(pupl$ploidy/2)*2==4,T,F)
  SAMPLEn=pupl$ploidy
  print(paste(SAMPLEpurity,SAMPLEn))
  # Estimating LogR deviation in diploid and gained regions (AUTOSOMAL)
  BB=read.table(paste0(tumourname,"_copynumber_extended.txt"),header=T,stringsAsFactors = F)
  
  BBdip=BB[which(BB$nMaj1_A==1 & BB$nMin1_A==1 & BB$frac1_A==1),]
  # correction for LogR values
  BBcorr=-mean(BBdip$LogR) #diploid only
  if (nrow(BBdip)<=1){
    print("likely WGD sample")
    BBdip=BB[which(BB$nMaj1_A==2 & BB$nMin1_A==2 & BB$frac1_A==1),]
    cnloh=BB[which(BB$nMaj1_A==2 & BB$nMin1_A==0 & BB$frac1_A==1),]
    if (nrow(cnloh)>0){
      BBcorr=-mean(cnloh$LogR)
    } else if (nrow(cnloh)==0){
      print("CRUDE estimation of BBcorr based on assumption of 2 copies vs ploidy")
      BBcorr=-log2(2/SAMPLEn)
    }
  }
  BBg1=BB[which(BB$nMaj1_A==2 & BB$nMin1_A==1 & BB$frac1_A==1),]
  BBg2=BB[which(BB$nMaj1_A==3 & BB$nMin1_A==1 & BB$frac1_A==1),]
  BBg3=BB[which(BB$nMaj1_A==4 & BB$nMin1_A==1 & BB$frac1_A==1),]
  BBg4=BB[which(BB$nMaj1_A==3 & BB$nMin1_A==2 & BB$frac1_A==1),] # likely observed in WGD samples
  
  # get max gain N:
  BBcomb=rbind(BBdip,BBg1,BBg2,BBg3,BBg4)
  maxNMaj=max(BBcomb$nMaj1_A)
  
  # SD for LogR values - diploid and gain regions
  BBsd=c(sd(BBdip$LogR),sd(BBg1$LogR),sd(BBg2$LogR),sd(BBg3$LogR))
  #BBsd_mean=mean(BBsd,na.rm=T)
  BBsd_max=max(BBsd, na.rm=T)
  BBsd_max=max(BBsd_max,0.05) # accept a minimum of 5% sd in LogR variation
  
  # BB LOH - estimating sd for LOH/loss events
  BBloh=BB[which(BB$nMaj1_A==1 & BB$nMin1_A==0 & BB$frac1_A==1),]
  if (nrow(BBloh)<=1){ #sd would be NA
    print("likely WGD sample or no clonal LOH event or just one single LOH event observed")
    BBloh=BB[which(BB$nMin1_A==0 & BB$frac1_A==1),] # all LOH events with varying nMaj1_A including 2:0 events
  }
  
  # expected ChrX logR values
  explogrgainX=function(x){log2((SAMPLEpurity*x+(1-SAMPLEpurity)*1)/1)}
  explogrGain=sapply(2:10000,explogrgainX) # up to 10000 copies!
  
  explogrLoss=max(log2(0+(1-SAMPLEpurity)*1),log2(0.01)) # if purity ~ 1, then purity of 0.99 is assumed for a realistic explogR estimate
  
  # assign CN
  SEG=data.frame()
  for (j in 1:nrow(SAMPLEsegs)){
    seg=SAMPLEsegs[j,]
    seg$type=ifelse(seg$mean<0,"loss","gain")
    
    # is segment different from zero?
    seg$mean=seg$mean+BBcorr
    
    if (seg$type=="gain"){
      seg$CNA=ifelse(seg$mean>(0+1.96*BBsd_max),"yes","no")
    } else {
      seg$CNA=ifelse(seg$mean<(0-1.96*BBsd_max),"yes","no")
    }
    # copy number
    if (seg$CNA=="yes"){
      if (seg$type=="gain"){
        rank=which(sort(c(explogrGain,seg$mean))==seg$mean) # rank of observed logR mean for segment among the expected logR values
        seg$CN=rank+1
        # clonality test
        if (rank==1){
          seg$clonal=ifelse(round(explogrGain[rank]-seg$mean,digits=2)<=round((BBsd_max/explogrGain[rank]),digits=2),"yes","no") # CV
          
        } else if (rank>=5){ # STOPS calling 'subclonal' events when copy number is >=5
          if (abs(seg$mean-explogrGain[rank-1])<abs(seg$mean-explogrGain[rank])){
            
            seg$clonal="yes"
            seg$CN=seg$CN-1
          } else {
            seg$clonal="yes"
          }
        } else {
          if (abs(seg$mean-explogrGain[rank-1])<abs(seg$mean-explogrGain[rank])){
            
            seg$clonal=ifelse(round(seg$mean-explogrGain[rank-1],digits = 2)<=round((BBsd_max/explogrGain[rank-1]),digits=2),"yes","no")
            if (seg$clonal=="yes"){
              seg$CN=seg$CN-1
            }
          }
          else {
            seg$clonal=ifelse(round(explogrGain[rank]-seg$mean,digits=2)<round((BBsd_max/explogrGain[rank]),digits = 2),"yes","no")
          }
        }
      } else if (seg$type=="loss"){
        seg$CN=0
        if (nrow(BBloh)>1){
          seg$clonal=ifelse(round(abs(explogrLoss-seg$mean),digits=2)<round(abs(sd(BBloh$LogR)/explogrLoss),digits=2),"yes","no")
        } else if (nrow(BBloh)<=1){ #sd would be NA
          seg$clonal=ifelse(round(abs(explogrLoss-seg$mean),digits=2)<round(abs(BBsd_max/explogrLoss),digits=2),"yes","no")
        }
      }
    }
    else {
      seg$CN=1
      seg$clonal=NA
      print(paste("no CNA for segment",j))
    }
    if (seg$arm=="p" & seg$end.pos>x_centromere[1]-1e6 & seg$CNA=="yes" & seg$end.pos<seg$start.pos+1e6){
      print("segment is p-arm centromere noise")
      print(seg)
    } else if (seg$arm=="q" & seg$end.pos<x_centromere[2]+1e6 & seg$CNA=="yes" & seg$end.pos<seg$start.pos+1e6){
      print("segment is q-arm centromere noise")
      print(seg)
    } else {
      SEG=rbind(SEG,seg)
    }
  }
  
  # CALCULATE CCF 
  CCF=data.frame()
  for (j in 1:nrow(SEG)){
    seg=SEG[j,]
    if (seg$CNA=="yes"){
      if (seg$type=="gain"){
        if (seg$clonal=="no"){
          seg$CCF=(2^seg$mean-(SAMPLEpurity*(seg$CN-1)+(1-SAMPLEpurity)*1))/SAMPLEpurity 
        } else {
          seg$CCF=1
        }
      } else if (seg$type=="loss"){
        if (seg$clonal=="no"){
          seg$CCF=(1-2^(seg$mean))/SAMPLEpurity # for Loss (assuming one chrX in all cells prior to Loss)
          if (seg$CCF>=0.95){
            seg$CCF=1
            seg$clonal="yes"
          }
        } else {
          seg$CCF=1
        }
      }
    } else {
      seg$CCF=1
    }
    CCF=rbind(CCF,seg)
  }
  
  # GENERATE FINAL OUTPUT
  SUBCLONES=data.frame()
  for (j in 1:nrow(CCF)){
    subclones=CCF[j,]
    if (subclones$CNA=="no"){
      subclones=data.frame(subclones,nMaj1=1,nMin1=0,frac1=1,nMaj2=0,nMin2=0,frac2=0)
    } else {
      if (subclones$type=="gain" & subclones$clonal=="yes"){
        subclones=data.frame(subclones,nMaj1=subclones$CN,nMin1=0,frac1=1,nMaj2=0,nMin2=0,frac2=0)  
      }
      else if (subclones$type=="gain" & subclones$clonal=="no"){
        if(subclones$CCF>0.5){ # switch nMaj/nMin so that the first nMaj/nMin represent the MAJOR CLONE
          subclones=data.frame(subclones,nMaj1=subclones$CN,nMin1=0,frac1=subclones$CCF,nMaj2=subclones$CN-1,nMin2=0,frac2=1-subclones$CCF)
        } else {
          subclones=data.frame(subclones,nMaj1=subclones$CN-1,nMin1=0,frac1=1-subclones$CCF,nMaj2=subclones$CN,nMin2=0,frac2=subclones$CCF)  
        }
      }
      else if (subclones$type=="loss" & subclones$clonal=="yes"){
        subclones=data.frame(subclones,nMaj1=subclones$CN,nMin1=0,frac1=1,nMaj2=0,nMin2=0,frac2=0) # very unlikely scenario; no sequencing reads should be present!
      }
      else if (subclones$type=="loss" & subclones$clonal=="no"){
        if(subclones$CCF>0.5){ # switch nMaj/nMin so that the first nMaj/nMin represent the MAJOR CLONE
          subclones=data.frame(subclones,nMaj1=subclones$CN,nMin1=0,frac1=subclones$CCF,nMaj2=1,nMin2=0,frac2=1-subclones$CCF)
        } else {
          subclones=data.frame(subclones,nMaj1=1,nMin1=0,frac1=1-subclones$CCF,nMaj2=subclones$CN,nMin2=0,frac2=subclones$CCF)
        }
      }
    }
    #print(j)
    SUBCLONES=rbind(SUBCLONES,subclones)
  }
  
  SUBCLONES$average=(SUBCLONES$nMaj1+SUBCLONES$nMin1)*SUBCLONES$frac1+(SUBCLONES$nMaj2+SUBCLONES$nMin2)*SUBCLONES$frac2
  
  SUBCLONESout=data.frame(SUBCLONES[,c("chrom","arm")],startpos=SUBCLONES$start.pos,endpos=SUBCLONES$end.pos,nSNPs=SUBCLONES$n.probes,
                          LogR=SUBCLONES$mean,SUBCLONES[,c("type","CNA","CN","clonal","nMaj1","nMin1","frac1","nMaj2","nMin2","frac2")],
                          subclonalCN=SUBCLONES$average,stringsAsFactors = F)
  SUBCLONESout$type[SUBCLONESout$type=="gain"]="+ve"
  SUBCLONESout$type[SUBCLONESout$type=="loss"]="-ve"
  
  # merge adjacent segments with same copy number  
  SUBCLONESout$rank=1:nrow(SUBCLONESout)
  SUBCLONESout=SUBCLONESout[order(SUBCLONESout$subclonalCN),]
  
  SPLIT=split(SUBCLONESout$rank, cumsum(c(1, diff(SUBCLONESout$rank) != 1))) # find consecutive segments with same subclonalCN
  outputDF=data.frame()
  for (j in 1:length(SPLIT)){
    if (length(SPLIT[[j]])>1){
      #print(length(SPLIT[[j]]))
      SUBsplit=SUBCLONESout[which(!is.na(match(SUBCLONESout$rank,SPLIT[[j]]))),]
      if (length(unique(SUBsplit$arm))==1){
        if (sd(SUBsplit$subclonalCN)<=0.01){
          mergedseg=SUBsplit[1,]
          mergedseg$endpos=SUBsplit[length(SPLIT[[j]]),"endpos"]
          mergedseg$nSNPs=sum(SUBsplit$nSNPs)
          mergedseg$LogR=weighted.mean(SUBsplit$LogR,SUBsplit$nSNPs)
          outputDF=rbind(outputDF,mergedseg)
        } else {
          outputDF=rbind(outputDF,SUBsplit)
          print("adjacent not same subclonalCN in SPLIT")
        }
      } else if (length(SPLIT[[j]])==2){
        outputDF=rbind(outputDF,SUBsplit)
      } else{
        # if (length(SUBsplit$arm=="p"))
        pseg=SUBsplit[SUBsplit$arm=="p",]
        if (nrow(pseg)>1){
          if (sd(pseg$subclonalCN)<=0.01){
            mergedseg=pseg[1,]
            mergedseg$endpos=pseg[nrow(pseg),"endpos"]
            mergedseg$nSNPs=sum(pseg$nSNPs)
            mergedseg$LogR=weighted.mean(pseg$LogR,pseg$nSNPs)
            outputDF=rbind(outputDF,mergedseg)
          } else {
            outputDF=rbind(outputDF,pseg)
            print("adjacent not same subclonalCN in pseg")
          }
        } else {outputDF=rbind(outputDF,pseg)}
        qseg=SUBsplit[SUBsplit$arm=="q",]
        if (nrow(qseg)>1){
          if (sd(qseg$subclonalCN)<=0.01){
            mergedseg=qseg[1,]
            mergedseg$endpos=qseg[nrow(qseg),"endpos"]
            mergedseg$nSNPs=sum(qseg$nSNPs)
            mergedseg$LogR=weighted.mean(qseg$LogR,qseg$nSNPs)
            outputDF=rbind(outputDF,mergedseg)
          } else {
            outputDF=rbind(outputDF,qseg)
            print("adjacent not same subclonalCN in qseg")
          }
        } else {outputDF=rbind(outputDF,qseg)}
      }
    }  else {
      SUBsplit=SUBCLONESout[which(SUBCLONESout$rank==SPLIT[[j]]),]
      outputDF=rbind(outputDF,SUBsplit)
    }
  }
  outputDF=outputDF[order(outputDF$startpos),]
  
  print(paste("Number of rows merged =",nrow(SUBCLONESout)-nrow(outputDF)))
  
  BBnew=BB[which(is.na(match(BB$chr,c("X","chrX")))),c(1:3,8:13)] # copynumber.txt columns to be populated with chrX calls
  
  outputDF_for_merge=data.frame(chr=outputDF$chrom,startpos=outputDF$startpos,endpos=outputDF$endpos,
                                nMaj1_A=outputDF$nMaj1,nMin1_A=outputDF$nMin1,frac1_A=outputDF$frac1,
                                nMaj2_A=outputDF$nMaj2,nMin2_A=outputDF$nMin2,frac2_A=outputDF$frac2,
                                stringsAsFactors = F)
  
  BBnew=rbind(BBnew,outputDF_for_merge)
  write.table(BBnew,paste0(tumourname,"_copynumber.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  
  BBnew_extended=BB[which(is.na(match(BB$chr,c("X","chrX")))),] # copynumber_extended.txt columns for chrX
  
  outputDF_for_merge_extended=data.frame(chr=outputDF$chrom,startpos=outputDF$startpos,endpos=outputDF$endpos,BAF=NA,pval=NA,LogR=outputDF$LogR,ntot=NA,
                                         nMaj1_A=outputDF$nMaj1,nMin1_A=outputDF$nMin1,frac1_A=outputDF$frac1,nMaj2_A=outputDF$nMaj2,nMin2_A=outputDF$nMin2,
                                         frac2_A=outputDF$frac2)
  BtoFsolutions=data.frame(matrix(nrow= nrow(outputDF),ncol = ncol(BB)-ncol(outputDF_for_merge_extended)))
  names(BtoFsolutions)=names(BB)[(ncol(outputDF_for_merge_extended)+1):ncol(BB)]
  
  BBnew_extended=rbind(BBnew_extended,cbind(outputDF_for_merge_extended,BtoFsolutions))
  write.table(BBnew_extended,paste0(tumourname,"_copynumber_extended.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  
  # PLOT
  outputDF$diff=outputDF$endpos-outputDF$startpos
  # print(outputDF)
  if (nrow(outputDF[which(outputDF$CNA=="yes"),])>0){
  PGAclonal=sum(outputDF[which(outputDF$clonal=="yes"),]$diff)/sum(outputDF[which(!is.na(outputDF$clonal)),]$diff)
  print(paste("chrX-based PGA.is.clonal =",PGAclonal))
  } else {
    print("no chrX CNA identified")
    PGAclonal = "NA"
  }
  
  plot_BB=ggplot()+geom_hline(yintercept = 0:ceiling(max(outputDF$subclonalCN)),linetype="longdash",col="grey",linewidth=0.2)+
    geom_rect(data=outputDF,aes(xmin=startpos,xmax=endpos,ymin=subclonalCN-0.02,ymax=subclonalCN+0.02))+
    geom_vline(xintercept = x_centromere,linetype="longdash",col="green")+
    #geom_hline(yintercept = nonpar,linetype="dotted",col="blue")+
    ylim(-0.2,ceiling(max(outputDF$subclonalCN))+0.2)+labs(x="ChrX coordinate (bp)",y="Average Ploidy")+
    theme(plot.title = element_text(hjust = 0.5,size=12),panel.background = element_blank())+
    ggtitle(paste0(tumourname," , Ploidy: ",round(SAMPLEn,digits = 3)," , Purity: ",round(SAMPLEpurity*100,digits = 0),
                   "%, chrX PGA.is.clonal: ",ifelse(PGAclonal=="NA","NA",paste0(round(PGAclonal*100,digits = 1),"%"))))
  
  # ANDROGEN RECEPTOR LOCUS
  if (AR){
    data.table::setDT(ar)
    data.table::setkey(ar,"startpos","endpos")
    data.table::setDT(outputDF)
    data.table::setkey(outputDF,"startpos","endpos")
    segAR=data.table::foverlaps(ar,outputDF,type="any",nomatch = 0)
    segAR$subclonalCN=(segAR$nMaj1+segAR$nMin1)*segAR$frac1+(segAR$nMaj2+segAR$nMin2)*segAR$frac2
    plot_BB=plot_BB+geom_rect(data=segAR,aes(xmin=startpos,xmax=endpos,ymin=subclonalCN-0.02,ymax=subclonalCN+0.02),fill="red")
  }
  
  pdf(paste0(tumourname,"_chrX_average_ploidy.pdf"))
  print(plot_BB)
  dev.off()
  
  # update outputDF (chrX-only copynumber output file)
  outputDF=outputDF[,c(1:6,11:17)]
  write.table(outputDF,paste0(tumourname,"_chrX_copynumber.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  
  # Update the genomewide Battenberg plots
  # goodness from rho_psi file (i.e. column named 'distance')
  goodness=read.table(paste0(tumourname,"_rho_and_psi.txt"),header=T,stringsAsFactors = F,sep="\t")
  goodness=goodness[which(goodness$is.best=="TRUE"),"distance"]
  # rho and ploidy from purity_ploidy file
  rho_psi=read.table(paste0(tumourname,"_purity_ploidy.txt"),header=T,stringsAsFactors = F,sep="\t")
  # update for BB3 - replace cellularity with purity
  # rho=rho_psi$cellularity
  rho=rho_psi$purity
  ploidy=rho_psi$ploidy
  # Need BAFsegment file
  BAFvals=as.data.frame(Battenberg:::read_bafsegmented(paste0(tumourname,".BAFsegmented.txt")))
  print("BAFvals")
  
  # replacing constant value of 90000 with chrX_BAFvals_length as a sample-specific way of counting the typical no. of het SNPs expected based on chrX length (chr 7 and 8 average hetSNP count) 
  # option 1 (may not always work if chr7 or chr8 have any kind of LOH in a pure or high-purity sample)
  #chrX_BAFvals_length=round((nrow(BAFvals[which(!is.na(match(BAFvals$Chromosome,c(7,"chr7")))),])+nrow(BAFvals[which(!is.na(match(BAFvals$Chromosome,c(8,"chr8")))),]))/2,0)
  # option 2 (based on the proportion of genome covered by chrX (i.e. 156e6/3e9 = 5%) and the number of hetSNPs in a sample-specific manner)
  chrX_BAFvals_length = round(nrow(BAFvals)*0.05,0)
  print(paste("chrX BAFvals length =",chrX_BAFvals_length))
  
  
  BAFvals=rbind(BAFvals[which(is.na(match(BAFvals$Chromosome,c("X","chrX")))),],
                data.frame(Chromosome="X",Position=sort(sample(1:155e6,chrX_BAFvals_length,replace=F)), # 155e6: approximate length of chrX
                           BAF=sample(c(0,1),chrX_BAFvals_length,replace=T),BAFphased=1,BAFseg=1)) 
  
  Battenberg:::plot.gw.subclonal.cn(subclones=BBnew, 
                                    BAFvals=BAFvals, 
                                    rho=rho, 
                                    ploidy=ploidy, 
                                    goodness=goodness, 
                                    output.gw.figures.prefix=paste(tumourname,"_BattenbergProfile", sep=""), 
                                    chr.names=chrom_names, 
                                    tumourname=tumourname)
}
