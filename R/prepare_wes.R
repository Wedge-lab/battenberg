
#' Add proxy SNPs for hetSNPs in the normal sample, impute their allelecounts and 
#' creating sample-specific alleles files matching the updated allelecounts files
#'
#'
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param tumourname Identifier to be used for tumour output files
#' @param normalname Identifier to be used for normal output files
#' @param g1000prefix_al Prefix path to the 1000 Genomes alleles reference files
#' @param proxy_snps Prefix path to the proxy SNP reference files generated for WES Battenberg based on high linkage disequilibrium (LD R2>=0.8) in the 'entire' 1000 Genomes SNP dataset for both hg19 and hg38
#' @param R2threshold Minimum linkage disequilibrium R2 value threshold for choosing proxy SNPs (range accepted = 0.8-1 inclusive, default 1 i.e. complete LD)
#' @param min_normal_depth Minimum depth required in the normal for a SNP to be included
#' @param heterozygousFilter The cutoff where a SNP in the normalname sample will be considered as heterozygous (default 0.1)
#' @author naser.ansari-pour
#' @export

add_snp_proxies=function(chrom_names,tumourname,normalname,g1000prefix_al,proxy_snps,R2threshold=1,min_normal_depth=10,heterozygousFilter=0.1){
  
  proxy_summary=data.frame()
  
  for (chrom in chrom_names){
    
    print(chrom)
    
    ac=read.table(paste0(normalname,"_alleleFrequencies_chr",chrom,".txt"),header=F,stringsAsFactors=F)
    names(ac)=c("Chr","Position","A","C","G","T","Depth")
    
    # Chr Notation standardisation within function
    ac$Chr=gsub("chr","",ac$Chr)
    ac$Chr=as.numeric(ac$Chr)
    
    al=read.table(paste0(g1000prefix_al,chrom,".txt"),header=T,stringsAsFactors=F)
    dim(al)
    ac$baf=ac[cbind(1:nrow(ac),al$a1+2)]/(ac[cbind(1:nrow(ac),al$a1+2)]+ac[cbind(1:nrow(ac),al$a0+2)])
    het_ac=ac[!is.nan(ac$baf),]
    het_ac=het_ac[which(het_ac$Depth>=min_normal_depth),]
    het_ac=het_ac[which(het_ac$baf>heterozygousFilter & het_ac$baf<(1-heterozygousFilter)),]
    
    proxy_snp_file=read.table(paste0(proxy_snps,chrom,".txt"),header=T,stringsAsFactors=F)
    proxy_snp_file=proxy_snp_file[which(proxy_snp_file$R2>=R2threshold),]
    
    # add proxy SNPs to the normal heterozygous SNPs and update the allelecounts by imputation
    
    het_ac_proxy=merge(het_ac,proxy_snp_file,by=c("Chr","Position"))
    
    # Find the rows with the maximum Depth for each unique Proxy SNP and remove other duplicated rows (rule: one additional row of imputed allelecounts per added Proxy SNP)
    if (nrow(het_ac_proxy)>0){
    het_ac_proxy$maxID=paste(het_ac_proxy$Proxy_Position,het_ac_proxy$Depth,sep="_")
    het_ac_proxy_unique=aggregate(Depth ~ Proxy_Position, data = het_ac_proxy, FUN = function(x){max(x)})
    het_ac_proxy_unique$ID=paste(het_ac_proxy_unique$Proxy_Position,het_ac_proxy_unique$Depth,sep="_")
    het_ac_proxy=het_ac_proxy[match(het_ac_proxy_unique$ID,het_ac_proxy$maxID),]
    print(paste0("Proportion of het SNPs added by proxy = ",round(nrow(het_ac_proxy)/nrow(het_ac),2)*100,"%"))
    proxy_summary=rbind(proxy_summary,data.frame(chrom=chrom,hetN=nrow(het_ac),proxyN=nrow(het_ac_proxy),R2threshold=R2threshold,hetsnp_increase=paste0(round(nrow(het_ac_proxy)/nrow(het_ac),2)*100,"%")))
    
    # output for normal sample
    out_ac=rbind(ac,data.frame(Chr=het_ac_proxy$Proxy_Chr,Position=het_ac_proxy$Proxy_Position,het_ac_proxy[,c("A","C","G","T","Depth","baf")]))
    out_ac=out_ac[order(out_ac$Position),]
    out_ac$baf=NULL
    
    
    
    cat(paste("#CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth",sep="\t"),sep="\n",file=paste0(normalname,"_alleleFrequencies_chr",chrom,".txt"))
    write.table(out_ac, paste0(normalname,"_alleleFrequencies_chr",chrom,".txt"),col.names=F,row.names=F,quote=F,sep="\t",append=T) # write over existing ac files
    
    # Update tumour ac according to the additions in the normal sample
    
    tac=read.table(paste0(tumourname,"_alleleFrequencies_chr",chrom,".txt"),header=F,stringsAsFactors=F)
    names(tac)=c("Chr","Position","A","C","G","T","Depth")
    
    # Chr Notation standardisation within function
    tac$Chr=gsub("chr","",tac$Chr)
    tac$Chr=as.numeric(tac$Chr)
    
    # add het_ac_proxy Proxy SNPs to tac and add imputed allelecounts from tac
    
    tac_proxy=merge(tac,het_ac_proxy[,c("Chr","Position","Proxy_Chr","Proxy_Position")],by=c("Chr","Position"))
    
    # Find the rows with the maximum Depth for each unique Proxy SNP and remove other duplicated rows (rule: one additional row of imputed allelecounts per added Proxy SNP)
    tac_proxy$maxID=paste(tac_proxy$Proxy_Position,tac_proxy$Depth,sep="_")
    tac_proxy_unique=aggregate(Depth ~ Proxy_Position, data = tac_proxy, FUN = function(x){max(x)})
    tac_proxy_unique$ID=paste(tac_proxy_unique$Proxy_Position,tac_proxy_unique$Depth,sep="_")
    tac_proxy=tac_proxy[match(tac_proxy_unique$ID,tac_proxy$maxID),]
    nrow(tac_proxy)
    
    if (nrow(tac_proxy)==nrow(het_ac_proxy)){
      out_tac=rbind(tac,data.frame(Chr=tac_proxy$Proxy_Chr,Position=tac_proxy$Proxy_Position,tac_proxy[,c("A","C","G","T","Depth")]))
      out_tac=out_tac[order(out_tac$Position),]
    }
    cat(paste("#CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth",sep="\t"),sep="\n",file=paste0(tumourname,"_alleleFrequencies_chr",chrom,".txt"))
    write.table(out_tac, paste0(tumourname,"_alleleFrequencies_chr",chrom,".txt"),col.names=F,row.names=F,quote=F,sep="\t",append=T) # write over existing ac files
    
    # create a temporary alleles file to match the updated/longer allelecount files with proxies for getBAFandLogR, etc.
    
    al_proxy=data.frame(position=het_ac_proxy$Proxy_Position,al[match(het_ac_proxy$Position,al$position),][,-1])
    dim(al_proxy)
    out_al=rbind(al,al_proxy)
    out_al=out_al[order(out_al$position),]
    
    write.table(out_al,paste0(tumourname,"_WES_1000genomes_alleles_",genomebuild,"_chr",chrom,".txt"),col.names=T,row.names=F,quote=F,sep="\t")
    } else {
      # define out_al with al as no 'proxy' SNP is there to be added but need to generate the alleles file for this chromosome with consistent naming i.e. *_WES_1000genomes_alleles_*
      # ac and tac stay the same so no output required
      out_al=al
      write.table(out_al,paste0(tumourname,"_WES_1000genomes_alleles_",genomebuild,"_chr",chrom,".txt"),col.names=T,row.names=F,quote=F,sep="\t")
    }
  }
  
  write.table(proxy_summary,paste0(tumourname,"_proxy_hetSNP_summary.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  
}


#' Prepare WES data for LogR estimation and hetSNP haplotype construction
#'
#' This function performs part of the Battenberg WES pipeline: Counting alleles, identifying hetSNPs in the normalname sample,
#' adding high LD proxy SNPs to observed hetSNPs, creating sample-specific alleles files, calculating BAF and LogR
#' and performing GC content correction.
#'
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param tumourbam Full path to the tumour BAM file
#' @param normalbam Full path to the normal BAM file
#' @param tumourname Identifier to be used for tumour output files
#' @param normalname Identifier to be used for normal output files
#' @param g1000prefix Prefix path to the 1000 Genomes SNP loci reference files specific to WES regions
#' @param g1000allelesprefix Prefix path to the 1000 Genomes alleles reference files
#' @param heterozygousFilter The cutoff where a SNP in the normalname sample will be considered as heterozygous (default 0.1)
#' @param proxy_snps Prefix path to the proxy SNP reference files generated for WES Battenberg based on high linkage disequilibrium (LD R2>=0.8) in the 'entire' 1000 Genomes SNP dataset for both hg19 and hg38
#' @param R2threshold Minimum linkage disequilibrium R2 value threshold for choosing proxy SNPs (range accepted = 0.8-1 inclusive, default 1 i.e. complete LD)
#' @param g1000allelesprefix_with_proxies Prefix path to the sample-specific alleles reference files to be generated in working directory that will include proxy SNPs for BAF & LogR estimation 
#' @param gccorrectprefix Prefix path to GC content reference data
#' @param repliccorrectprefix Prefix path to replication timing reference data (supply NULL if no replication timing correction is to be applied)
#' @param min_base_qual Minimum base quality required for a read to be counted
#' @param min_map_qual Minimum mapping quality required for a read to be counted
#' @param allelecounter_exe Path to the allele counter executable (can be found in $PATH)
#' @param min_normal_depth Minimum depth required in the normal for a SNP to be included
#' @param nthreads The number of parallel processes to run for allelecounting
#' @author naser.ansari-pour
#' @export

prepare_wes = function(chrom_names, tumourbam, normalbam, tumourname, normalname, g1000prefix, g1000allelesprefix, heterozygousFilter,
                       proxy_snps, R2threshold,g1000allelesprefix_with_proxies,gccorrectprefix,repliccorrectprefix, min_base_qual, 
                       min_map_qual, allelecounter_exe, min_normal_depth) {
  
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  
  # Obtain allele counts for 1000 Genomes WES locations for both tumour and normal - must NOT be skipped because add_snp_proxies updates them
  foreach::foreach(i=1:length(chrom_names),.export="getAlleleCounts") %dopar% {
    getAlleleCounts(bam.file=tumourbam,
                      output.file=paste(tumourname,"_alleleFrequencies_chr", chrom_names[i], ".txt", sep=""),
                      g1000.loci=paste(g1000prefix, chrom_names[i], ".txt", sep=""),
                      min.base.qual=min_base_qual,
                      min.map.qual=min_map_qual,
                      allelecounter.exe=allelecounter_exe)
      
    getAlleleCounts(bam.file=normalbam,
                        output.file=paste(normalname,"_alleleFrequencies_chr", chrom_names[i], ".txt",  sep=""),
                        g1000.loci=paste(g1000prefix, chrom_names[i], ".txt", sep=""),
                        min.base.qual=min_base_qual,
                        min.map.qual=min_map_qual,
                        allelecounter.exe=allelecounter_exe)
      }
  
  # identify hetSNPs in normal and add proxy SNPs
  add_snp_proxies(tumourname = tumourname,
                  normalname = normalname,
                  g1000prefix_al = g1000allelesprefix,
                  proxy_snps = proxy_snps,
                  chrom_names = chrom_names)
  
  # Obtain BAF and LogR from the update allele counts with proxy SNPs + sample-specific alleles files generated in the working directory
  getBAFsAndLogRs(tumourAlleleCountsFile.prefix=paste(tumourname,"_alleleFrequencies_chr", sep=""),
                  normalAlleleCountsFile.prefix=paste(normalname,"_alleleFrequencies_chr", sep=""),
                  figuresFile.prefix=paste(tumourname, "_", sep=''),
                  BAFnormalFile=paste(tumourname,"_normalBAF.tab", sep=""),
                  BAFmutantFile=paste(tumourname,"_mutantBAF.tab", sep=""),
                  logRnormalFile=paste(tumourname,"_normalLogR.tab", sep=""),
                  logRmutantFile=paste(tumourname,"_mutantLogR.tab", sep=""),
                  combinedAlleleCountsFile=paste(tumourname,"_alleleCounts.tab", sep=""),
                  chr_names=chrom_names,
                  g1000file.prefix=g1000allelesprefix_with_proxies,
                  minCounts=min_normal_depth,
                  samplename=tumourname)
  
  # Perform GC correction - identical to WGS run
  gc.correct.wgs(Tumour_LogR_file=paste(tumourname,"_mutantLogR.tab", sep=""),
                 outfile=paste(tumourname,"_mutantLogR_gcCorrected.tab", sep=""),
                 correlations_outfile=paste(tumourname, "_GCwindowCorrelations.txt", sep=""),
                 gc_content_file_prefix=gccorrectprefix,
                 replic_timing_file_prefix=repliccorrectprefix,
                 chrom_names=chrom_names)
}
