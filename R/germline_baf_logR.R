germline_baf_logR <- function(GERMLINENAME,G1000PREFIX,CHROM_NAMES){
  chrom_names=CHROM_NAMES
  #read heterozygous SNPs per chromosome for alleleCounter files & 1000G allele files####
  AC=list() # alleleCounts
  AL=list() # 1000G alleles
  MaC=list() # matched alleleCounts
  OHET=list() # HET SNP data
  for (chr in chrom_names){
    # read in alleleCounter output for each chromosome
    ac=read.table(paste0(GERMLINENAME,"_alleleFrequencies_chr",chr,".txt"),stringsAsFactors = F)
    ac=ac[order(ac$V2),]
    AC[[chr]]=ac
    print(length(AC))
    # match allele counts with respective SNP alleles
    
    al=read.table(paste0(G1000PREFIX,chr,".txt"),header=T,stringsAsFactors = F)
    AL[[chr]]=al
    print(length(AL))
    #etc
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
    #extract rows with 0.1=<baf=<0.9 - final variable 'oo'
    ohet=o[which(o$baf>=0.10 & o$baf<=0.90 & o$depth>10),]
    # this output will replace 'd'
    #d=read.table("~/Documents/BDI/Oxford-Celgene_MM/D1_12_chr1_heterozygousMutBAFs_haplotyped.txt",header=T,stringsAsFactors = F)
    ohet$Position2=c(ohet$Position[2:nrow(ohet)],2*ohet$Position[nrow(ohet)]-ohet$Position[nrow(ohet)-1])
    ohet$Position_dist=ohet$Position2-ohet$Position
    ohet$Position_dist_percent=ohet$Position_dist/max(ohet$Position_dist)
    OHET[[chr]]=ohet
    print(paste("chromosome",chr,"file read"))
  }
  # CREATE mutantBAF and mutantLogR *.tab files #
  germline=GERMLINENAME
  MAC=data.frame()
  for (chr in chrom_names){
    MaC_CHR=data.frame(chr=chr,MaC[[chr]])
    MAC=rbind(MAC,MaC_CHR)
    print(chr)
  }
  names(MAC)=c("chr","position","a0","a1","ref","alt","coverage","baf")
  print(head(MAC))
  print(dim(MAC))
  MAC$logr=log2(MAC$coverage/mean(MAC$coverage))
  MACC=MAC[which(!is.na(MAC$baf)),]
  print(nrow(MAC)-nrow(MACC))
  
  BAF=data.frame(Chromosome=MACC$chr,Position=MACC$pos,germline=MACC$baf)
  names(BAF)[names(BAF) == "germline"] <- germline
  BAF=BAF[order(BAF$Chromosome,BAF$Position),]
  BAF$Chromosome[BAF$Chromosome==23]="X"
  write.table(BAF,paste0(germline,"_mutantBAF.tab"),col.names=T,row.names=F,quote=F,sep="\t")
  rm(BAF)
  
  LogR=data.frame(Chromosome=MACC$chr,Position=MACC$pos,germline=MACC$logr)
  names(LogR)[names(LogR) == "germline"] <- germline
  LogR=LogR[order(LogR$Chromosome,LogR$Position),]
  LogR$Chromosome[LogR$Chromosome==23]="X" # revert back from 23 to X for Chromosome number
  write.table(LogR,paste0(germline,"_mutantLogR.tab"),col.names=T,row.names=F,quote=F,sep="\t")
  
  rm(MAC)
  rm(MaC)
  rm(MACC)
  GL_OHET <<- OHET
  GL_AL <<- AL
  GL_AC <<- AC
  GL_LogR <<- LogR 
  print("STEP 1 - BAF and LogR - completed")
}
