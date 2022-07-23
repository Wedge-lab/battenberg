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
    } else {
      
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
      }
    }
  print("STEP 1 - BAF and LogR - completed")
}

