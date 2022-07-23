#' Obtain BAF and LogR for tumour-only mode allele counts
#'
#' Function to generate BAF and LogR files based on allele counts of the Cell line.
#' It also generates the input data required by the following 'cell_line_reconstruct_normal' function.
#' @param tumourname The tumour name used for Battenberg (i.e. the cell line BAM file name without the '.bam' extension).
#' @param g1000alleles.prefix Prefix to where 1000 Genomes allele files can be found.
#' @param chrom_names A vector with allowed chromosome names.
#' @param snv_rho Estimated purity or aberrant cell fraction of the tumour sample based on SNV VAF-based approach 
#' @author Naser Ansari-Pour (WIMM, Oxford)
#' @export

tumour_only_baf_logR = function(tumourname,g1000alleles.prefix,chrom_names,snv_rho){
  # read heterozygous SNPs per chromosome for alleleCounter files & 1000G allele files
  AC=list() # alleleCounts
  AL=list() # 1000G alleles
  MaC=list() # matched alleleCounts
  OHET=list() # HET SNP data
  for (chr in chrom_names){
    # read in alleleCounter output for each chromosome
    ac=read.table(paste0(tumourname,"_alleleFrequencies_chr",chr,".txt"),stringsAsFactors = F)
    ac=ac[order(ac$V2),]
    AC[[chr]]=ac
    print(length(AC))
    # match allele counts with respective SNP alleles
    al=read.table(paste0(g1000alleles.prefix,chr,".txt"),header=T,stringsAsFactors = F)
    AL[[chr]]=al
    print(length(AL))
    ref=al$a0
    ref_df=data.frame(pos=1:nrow(al),ref=ref+2)
    REF=ac[cbind(ref_df$pos,ref_df$ref)]
    alt=al$a1
    alt_df=data.frame(pos=1:nrow(al),alt=alt+2)
    ALT=ac[cbind(alt_df$pos,alt_df$alt)]
    mac=data.frame(ref=REF,alt=ALT)
    mac$depth=as.numeric(mac$ref)+as.numeric(mac$alt)
    mac$baf=as.numeric(mac$alt)/as.numeric(mac$depth)
    o=cbind(al,mac)
    names(o)=c("Position","a0","a1","ref","alt","depth","baf")
    MaC[[chr]]=o
  }
  # CREATE mutantBAF and mutantLogR *.tab files 
  MAC=data.frame()
  for (chr in chrom_names){
    MaC_CHR=data.frame(chr=chr,MaC[[chr]])
    MAC=rbind(MAC,MaC_CHR)
    print(chr)
  }
  names(MAC)=c("chr","position","a0","a1","ref","alt","coverage","baf")
  print(head(MAC))
  print(dim(MAC))
  MACC=MAC[which(!is.na(MAC$baf)),]
  MACC$logr=log2(MACC$coverage/mean(MACC$coverage,na.rm=TRUE)) # in case of coverage == NA due to non-matching alleles or presence of indels in loci file
  
  BAF=data.frame(Chromosome=MACC$chr,Position=MACC$pos,tumourname=MACC$baf)
  names(BAF)[names(BAF) == "tumourname"] <- tumourname
  BAF=BAF[order(BAF$Chromosome,BAF$Position),]
  BAF$Chromosome[BAF$Chromosome==23]="X" # revert back from 23 to X for Chromosome name
  write.table(BAF,paste0(tumourname,"_mutantBAF.tab"),col.names=T,row.names=F,quote=F,sep="\t")
  rm(BAF)
  
  LogR=data.frame(Chromosome=MACC$chr,Position=MACC$pos,tumourname=MACC$logr)
  names(LogR)[names(LogR) == "tumourname"] <- tumourname
  LogR=LogR[order(LogR$Chromosome,LogR$Position),]
  LogR$Chromosome[LogR$Chromosome==23]="X" # revert back from 23 to X for Chromosome name
  write.table(LogR,paste0(tumourname,"_mutantLogR.tab"),col.names=T,row.names=F,quote=F,sep="\t")
  
  baf_threshold=2/mean(MACC$coverage) # 2 = min no. of reads supporting the SNP (alt>=2)
  
  # output for tumour_only_generate_normal
  if (!is.na(snv_rho)){
    if (snv_rho<0.95){
    # get heterozygous filter for run_haplotyping_tumour_only
    heterozygousFilter <<- baf_threshold
    } else if (snv_rho>0.95 & snv_rho<=1){
      
      # extract hetSNPs by taking into account alt>=2, minCount=10 and baf range depending on average depth of loci with coverage > 0
      # (i.e. regions with mapped reads and not whole-genome average of the sample which includes all coverage==0)
      
      for (chr in chrom_names){
          ohet=MACC[which(MACC$chr==chr & MACC$coverage>=10 & MACC$baf>=baf_threshold & MACC$baf<=(1-baf_threshold)),]
          names(ohet)[names(ohet) == "position"]="Position"
          ohet$Position2=c(ohet$Position[2:nrow(ohet)],2*ohet$Position[nrow(ohet)]-ohet$Position[nrow(ohet)-1])
          ohet$Position_dist=ohet$Position2-ohet$Position
          ohet$Position_dist_percent=ohet$Position_dist/max(ohet$Position_dist)
          OHET[[chr]]=ohet
          print(paste("chromosome",chr,"hetSNP output read"))
        }
        
        rm(MAC)
        rm(MaC)
        rm(MACC)
        TUM_OHET <<- OHET
        TUM_AL <<- AL
        TUM_AC <<- AC
        TUM_LogR <<- LogR
        heterozygousFilter <<- baf_threshold
      }
    }
  print("STEP 1 - BAF and LogR - completed")
}

#' Prepare WGS tumour-only data for haplotype construction
#' 
#' This function performs part of the Battenberg WGS pipeline: Counting alleles, generating BAF and logR, 
#' reconstructing normal-pair allele counts for high-purity (rho > 0.95) tumours and performing GC content correction.
#' 
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param tumourbam Full path to the tumour BAM file 
#' @param tumourname Identifier to be used for tumour output files (i.e. the cell line BAM file name without the '.bam' extension).
#' @param g1000lociprefix Prefix path to the 1000 Genomes loci reference files
#' @param g1000allelesprefix Prefix path to the 1000 Genomes SNP allele reference files
#' @param snv_rho Estimated purity or aberrant cell fraction of the tumour sample based on SNV VAF-based approach 
#' @param gamma_ivd The PCF gamma value for segmentation of 1000G hetSNP IVD values (Default 1e5).
#' @param kmin_ivd The min number of SNPs to support a segment in PCF of 1000G hetSNP IVD values (Default 50)
#' @param centromere_noise_seg_size The maximum size of PCF segment to be removed as noise when it overlaps with the centromere due to the noisy nature of data (Default 1e6)
#' @param centromere_dist The minimum distance from the centromere to ignore in analysis due to the noisy nature of data in the vicinity of centromeres (Default 5e5)
#' @param min_het_dist The minimum distance for detecting higher resolution inter-hetSNP regions with potential LOH while accounting for inherent homozygote stretches (Default 1e5)
#' @param gamma_logr The PCF gamma value for confirming LOH within each inter-hetSNP candidate segment (Default 100)
#' @param length_adjacent The length of adjacent regions either side of a candidate inter-hetSNP LOH region to be plotted (Default 5e4)
#' @param gccorrectprefix Prefix path to GC content reference data
#' @param repliccorrectprefix Prefix path to replication timing reference data (supply NULL if no replication timing correction is to be applied)
#' @param min_base_qual Minimum base quality required for a read to be counted
#' @param min_map_qual Minimum mapping quality required for a read to be counted
#' @param allelecounter_exe Path to the allele counter executable (can be found in $PATH)
#' @param min_normal_depth Minimum depth required in the normal for a SNP to be included
#' @param skip_allele_counting Flag, set to TRUE if allele counting is already complete (files are expected in the working directory on disk)
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
prepare_wgs_tumour_only = function(chrom_names, chrom_coord, tumourbam, tumourname, g1000lociprefix, g1000allelesprefix, snv_rho, gamma_ivd=1e5, kmin_ivd=50, centromere_noise_seg_size=1e6, 
                                 centromere_dist=5e5, min_het_dist=1e5, gamma_logr=100, length_adjacent=5e4, gccorrectprefix,repliccorrectprefix, min_base_qual, min_map_qual, 
                                 allelecounter_exe, min_normal_depth, skip_allele_counting) {
  
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  
  if (!skip_allele_counting) {
    # Obtain allele counts for 1000 Genomes locations for the cell line
    foreach::foreach(i=1:length(chrom_names)) %dopar% {
      getAlleleCounts(bam.file=tumourbam,
                      output.file=paste(tumourname,"_alleleFrequencies_chr", i, ".txt", sep=""),
                      g1000.loci=paste(g1000lociprefix, i, ".txt", sep=""),
                      min.base.qual=min_base_qual,
                      min.map.qual=min_map_qual,
                      allelecounter.exe=allelecounter_exe)
    }
  }
  
  # Standardise Chr notation (removes 'chr' string if present; essential for cell_line_baf_logR)
  
  standardiseChrNotation(tumourname=tumourname,
                         normalname=NULL) 
  
  # Obtain BAF and LogR from the raw allele counts of the cell line
  tumour_only_baf_logR(tumourname=tumourname,
                     g1000alleles.prefix=g1000allelesprefix,
                     chrom_names=chrom_names,
                     snv_rho=snv_rho
  )
  
  MIN_RHO <<- snv_rho-0.01
  MAX_RHO <<- snv_rho+0.01
  
  if (snv_rho<=0.95){
    print("STEP 2 completed - Min and Max rho updated")
  } else if (snv_rho>0.95 & snv_rho<=1){
  # Reconstruct normal-pair allele count files for high-purity tumours
  
  foreach::foreach(i=1:length(chrom_names),.export=c("tumour_only_reconstruct_normal","TUM_OHET","TUM_AL","TUM_AC","TUM_LogR"),.packages=c("copynumber","ggplot2","grid")) %dopar% {
    
    tumour_only_reconstruct_normal(TUMOURNAME=tumourname,
                                 NORMALNAME=paste0(tumourname,"_normal"),
                                 chrom_coord=chrom_coord,
                                 chrom=i,
                                 TUM_OHET=TUM_OHET,
                                 TUM_AL=TUM_AL,
                                 TUM_AC=TUM_AC,
                                 TUM_LogR=TUM_LogR,
                                 GAMMA_IVD=gamma_ivd,
                                 KMIN_IVD=kmin_ivd,
                                 CENTROMERE_NOISE_SEG_SIZE=centromere_noise_seg_size,
                                 CENTROMERE_DIST=centromere_dist,
                                 MIN_HET_DIST=min_het_dist,
                                 GAMMA_LOGR=gamma_logr,
                                 LENGTH_ADJACENT=length_adjacent)
  }
  
  if (length(list.files(pattern="normal_alleleFrequencies"))==length(chrom_names)){
    print("STEP 2 - Normal allelecounts reconstruction - completed")
  } else { 
    stop("Missing 'normal' allelecount files - all chromosomes NOT reconstructed")
  }
  }
  # Perform GC correction
  gc.correct.wgs(Tumour_LogR_file=paste(tumourname,"_mutantLogR.tab", sep=""),
                 outfile=paste(tumourname,"_mutantLogR_gcCorrected.tab", sep=""),
                 correlations_outfile=paste(tumourname, "_GCwindowCorrelations.txt", sep=""),
                 gc_content_file_prefix=gccorrectprefix,
                 replic_timing_file_prefix=repliccorrectprefix,
                 chrom_names=chrom_names)
  print("GC and replication correction completed")
}
