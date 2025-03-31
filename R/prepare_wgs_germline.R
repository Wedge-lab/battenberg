#' Chromosome notation standardisation (removing 'chr' string from chromosome names - mainly an issue in hg38 BAMs)
#'
#' @param GERMLINENAME The germline identifier, this is used as a prefix for the allele count files. If allele counts are supplied separately, they are expected to have this identifier as prefix.
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
standardiseChrNotation_germline = function(GERMLINENAME) {
gAF=capture.output(cat('bash -c \'sed -i \'s/chr//g\' ', GERMLINENAME,'_alleleFrequencies_chr*.txt\'',sep = ""))
system(gAF)
}

#' Obtain BAF and LogR from the Germline allele counts
#'
#' Function to generate BAF and LogR files based on allele counts of the Germline.
#' It also generates the input data required by the following 'germline_reconstruct_normal' function.
#' @param GERMLINENAME The germline name used for Battenberg (i.e. the Germline BAM file name without the '.bam' extension).
#' @param g1000alleles.prefix Prefix to where 1000 Genomes allele files can be found.
#' @param chrom_names A vector with allowed chromosome names.
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export

germline_baf_logR = function(GERMLINENAME,g1000alleles.prefix,chrom_names){
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
    
    al=read.table(paste0(g1000alleles.prefix,chr,".txt"),header=T,stringsAsFactors = F)
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
    #extract rows with 0.1=<baf=<0.9
    ohet=o[which(o$baf>=0.10 & o$baf<=0.90 & o$depth>10),]
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
  # MAC$logr=log2(MAC$coverage/mean(MAC$coverage))
  MAC$logr=log2(MAC$coverage/mean(MAC$coverage,na.rm=TRUE)) # in case of coverage == NA due to non-matching alleles or presence of indels in loci file
  MACC=MAC[which(!is.na(MAC$baf)),]
  print(nrow(MAC)-nrow(MACC))
  
  BAF=data.frame(Chromosome=MACC$chr,Position=MACC$pos,germline=MACC$baf)
  names(BAF)[names(BAF) == "germline"] <- germline
  BAF=BAF[order(BAF$Chromosome,BAF$Position),]
  BAF$Chromosome[BAF$Chromosome==23]="X" # revert back from 23 to X for Chromosome number
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

#' Reconstruct normal-pair allele count files for Germlines
#'
#' Function to generate normal-pair allele count files based on IVD-PCF and inter-hetSNP logR-based LOH detection (IVD: Inter-Variant Distance, het: heterozygote) 
#' This method reconstructs the normal-pair counts by using the allele counts of the Germline as template
#' It fills the detected LOH regions with evenly-distributed hetSNPs with the density estimated based on each chromosome in each germline sample
#' It essentially informs Battenberg of the location of hetSNPs across the genome in the germline sample
#' @param GERMLINENAME The germline name used for Battenberg (i.e. the germline BAM file name without the '.bam' extension)
#' @param NORMALNAME The normal name used for naming the generated normal-pair allele counts files
#' @param chrom_coord Full path to the file with chromosome coordinates including start, end and left/right centromere positions
#' @param chrom Chromosome number for which normal-pair will be reconstructed
#' @param GL_OHET List of observed heterozygous SNPs across all chromosomes generated within the germline_baf_logR function
#' @param GL_AL List of alleles at SNPs across all chromosomes generated within the germline_baf_logR function
#' @param GL_AC List of allele counts at SNPs across all chromosomes generated within the germline_baf_logR function
#' @param GL_LogR Dataframe of genomewide LogR values for SNPs across all chromosomes generated within the germline_baf_logR function
#' @param GAMMA_IVD The PCF gamma value for segmentation of 1000G hetSNP IVD values (Default 1e5)
#' @param KMIN_IVD The min number of SNPs to support a segment in PCF of 1000G hetSNP IVD values (Default 50)
#' @param CENTROMERE_DIST The minimum distance from the centromere to ignore in analysis due to the noisy nature of data in the vicinity of centromeres (Default 5e5)
#' @param MIN_HET_DIST The minimum distance for detecting higher resolution inter-hetSNP regions with potential LOH while accounting for inherent homozygote stretches (Default 1e5)
#' @param GAMMA_LOGR The PCF gamma value for confirming LOH within each inter-hetSNP candidate segment (Default 100)
#' @param LENGTH_ADJACENT The length of adjacent regions either side of a candidate inter-hetSNP LOH region to be plotted (Default 5e4)
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export

germline_reconstruct_normal = function(GERMLINENAME,NORMALNAME,chrom_coord,chrom,GL_OHET,GL_AL,GL_AC,GL_LogR,GAMMA_IVD,KMIN_IVD,CENTROMERE_NOISE_SEG_SIZE,CENTROMERE_DIST,MIN_HET_DIST,GAMMA_LOGR,LENGTH_ADJACENT){
  # IDENTIFY REGIONS OF LOH #
  colClasses=c(chr="numeric",start="numeric",cen.left.base="numeric",cen.right.base="numeric",end="numeric")
  chr_loc=read.table(chrom_coord,colClasses = colClasses,header=T,stringsAsFactors = F) # chrom_coord = full path to chromosome coordinates
  chr_loc$length=(chr_loc$cen.left.base-chr_loc$start)+(chr_loc$end-chr_loc$cen.right.base)
  #STEP 2.0: identify LOH by IVD-PCF
  LOH=list()
  PCF_folder = "PCF_plots"
  if(!file.exists(PCF_folder)){
    dir.create(PCF_folder)
  }
  i=chrom
  print(paste("chrom=",i))
  pcf_input=data.frame(chr=i,position=GL_OHET[[i]]$Position,IVD=(GL_OHET[[i]]$Position_dist_percent))
  pcf_input=pcf_input[which(pcf_input$position<chr_loc[i,"cen.left.base"]-CENTROMERE_DIST | pcf_input$position>chr_loc[i,"cen.right.base"]+CENTROMERE_DIST),]
  pcf_input=pcf_input[which(pcf_input$position>=chr_loc[i,"start"] & pcf_input$position<=chr_loc[i,"end"]),] # use only regions covered with gcCorrect LogR range 
  PCF=pcf(pcf_input,gamma=GAMMA_IVD,kmin = KMIN_IVD)
  pdf(paste0(PCF_folder,"/",GERMLINENAME,"_chr",i,"_PCF_plot.pdf"))
  plotChrom(pcf_input,PCF)
  dev.off()
  PCF$diff=PCF$end.pos-PCF$start.pos
  
  # Decide if there is any LOH based on PCF and chr_snp_density
  chr_snp_density=nrow(pcf_input)/(pcf_input$position[nrow(pcf_input)]-pcf_input$position[1]) # density of HET SNPs across the region covered by HET SNPs
  #CALCULATE min_normal_snp_density#
  # minimum normal density for SNPs (in bps) is 3 x 10^-4 with median of 7 x 10^-4
  ####
  min_normal_snp_density=0.0001
  loh_regions=PCF[which(round(PCF$mean,3)>0.001),] # LOH regions
  loh_regions=loh_regions[which(loh_regions$n.probes>1),] # only keep segments with minimum of 2 probes (SNPs) in PCF jump
  if (nrow(loh_regions)>0){
  if (mean(pcf_input$IVD)>0.01 & chr_snp_density<min_normal_snp_density){ #can change chr_snp_density from 0.00005 to 0.0001 as conservative measure - done
    #mean(pcf_input$IVD) or mean(PCF$mean) indicates presence of jumps in IVD
    loh_regions=loh_regions # LOH regions
    print(paste("full-length chromosomal loss at chr",i))
  } else if (sum(loh_regions$diff)>=((pcf_input$position[nrow(pcf_input)]-pcf_input$position[1]))*0.9 & chr_snp_density>min_normal_snp_density){
    # do PCF regions cover >=90% of the chromosome & is the chromosome snp density above the minimum
    loh_regions=0 # LOH regions
    print(paste("no PCF jumps at chr",i))
  } else {
    loh_regions=loh_regions # LOH regions
    print(paste("likely partial LOH(s) at chr",i))
  }
  } else {loh_regions=0}

  # loop to turn empty dataframe to 0 for loh_regions 
  #suppressWarnings(
  #  if (loh_regions[1]!=0){
  #   if (nrow(loh_regions)==0){
  #      loh_regions=0
  #    } else {print("dataframe non-empty")}
  #  } else {print("no LOH at all")})
  
  #filter regions for those next to the centromere and 'short'
  noise=NULL
  if (!is.null(nrow(loh_regions))){
    for (j in 1:nrow(loh_regions)){
      if (loh_regions$arm[j]=="p"){
        #if (loh_regions$end.pos[j]-chr_loc$cen.left.base[i]<1e5 & loh_regions$diff[j]<1e6){ #FOR EXCLUSION: max distance to centromere = 100kb , max length of short LOH region = 1Mb
        #  noise=append(noise,j)
        #}
        if (loh_regions$end.pos[j]>chr_loc$cen.left.base[i] & loh_regions$diff[j]<CENTROMERE_NOISE_SEG_SIZE){ #FOR EXCLUSION: segment is short IVD region (default<1Mb) and endpos is over the p-arm limit (ending point)
          noise=append(noise,j)
        }
        #if (loh_regions$end.pos[j]>chr_loc$cen.left.base[i] & loh_regions$diff[j]>CENTROMERE_NOISE_SEG_SIZE & !is.na(match(chrom,c(1,9,16)))){ # Chr 1,9,16 have large heterochromatin region next to centromere
        #  noise=append(noise,j)
        #}
      }
      if (loh_regions$arm[j]=="q"){
        #if (loh_regions$start.pos[j]-chr_loc$cen.right.base[i]<1e5 & loh_regions$diff[j]<1e6){ #FOR EXCLUSION: max distance to centromere = 100kb , max length of short LOH region = 1Mb
        #  noise=append(noise,j)
        #}
        if (loh_regions$start.pos[j]<chr_loc$cen.right.base[i] & loh_regions$diff[j]<CENTROMERE_NOISE_SEG_SIZE){ #FOR EXCLUSION: segment is short IVD region (default<1Mb) and startpos is below the q-arm limit (starting point)
          noise=append(noise,j)
        }
        if (loh_regions$start.pos[j]<(chr_loc$cen.right.base[i]+1e5) & loh_regions$diff[j]>CENTROMERE_NOISE_SEG_SIZE & !is.na(match(chrom,c(1,9,16)))){ # qARM of Chr 1,9,16 have large heterochromatin region next to centromere + 100kb tolerance for start of heterochromatin region
          noise=append(noise,j)
        }
      }
    }
  } else {print("no 'centromere noise' calculation")}
  if (!is.null(noise)){
    LOH_regions=loh_regions[-noise,]
  } else {LOH_regions=loh_regions}
  
  #remove LOH regions in the p arm of acrocentric chromosomes 13,14,15,21 and 22
  if (!is.na(match(i,c(13:15,21:22))) & !is.null(nrow(LOH_regions))){
    LOH_regions=LOH_regions[which(LOH_regions$arm!="p"),]
  }
  #
  if (is.null(dim(LOH_regions))){
    print(paste("no LOH detected in chr",i))
    LOH[[i]]=0
  } else if (dim(LOH_regions)[1]!=0 & dim(LOH_regions)[2]!=0) {
    print(paste("we have LOH for",sum(LOH_regions$diff),"bp in chr",i))
    LOH[[i]]=data.frame(chr=i,LOH_regions)
  } else if (dim(LOH_regions)[1]==0) {
    print(paste("no LOH regions remained after noise correction for chr",i))
    LOH[[i]]=0
  } else {print("unkown issue!")}
  print(paste("chrom=",i,"IVD-PCF finished"))
  #
  ##
  # STEP 2 - get higher resolution LOH regions
  ##
  #
  print(paste("chrom=",i))
  # use loop to find blocks with no LOH - while taking account of the centromere - RUN1
  ac=GL_AC[[i]]
  al=GL_AL[[i]]
  names(ac)=c("chr","position",1:4,"depth")
  chr_interval=c(chr_loc[i,"start"],chr_loc[i,"end"]) # use gcCorrect LogR range for chromosome interval
  if (!is.null(nrow(LOH[[i]]))){
    non_LOH=data.frame()## get all non_LOH regions ## 
    for (j in 1:(nrow(LOH[[i]])+1)){
      if (j == 1 & chr_interval[1]==LOH[[i]]$start.pos[j]){
        print("LOH from start of chromosome")
      } else if (j == 1 & chr_interval[1]<LOH[[i]]$start.pos[j]){
        non_loh=data.frame(start=chr_interval[1],end=LOH[[i]]$start.pos[j]-1)
      } else if (j>1 & j <= nrow(LOH[[i]]) & LOH[[i]]$arm[j]==LOH[[i]]$arm[j-1]){
        non_loh=data.frame(start=LOH[[i]]$end.pos[j-1]+1,end=LOH[[i]]$start.pos[j]-1)
      } else if (j>1 & j <= nrow(LOH[[i]]) & LOH[[i]]$arm[j]!=LOH[[i]]$arm[j-1]){
        non_loh=data.frame(start=c(LOH[[i]]$end.pos[j-1]+1,chr_loc[i,]$cen.right.base),end=c(chr_loc[i,]$cen.left.base,LOH[[i]]$start.pos[j]-1))
      } else{
        if ((LOH[[i]]$end.pos[j-1]+1)<chr_interval[2]){ # avoids going over the chromosome interval
          non_loh=data.frame(start=LOH[[i]]$end.pos[j-1]+1,end=chr_interval[2])
        } else{
          print("reached end of chromosome")
          rm(non_loh)
        }
      }
      print(j)
      if (exists("non_loh")){
        non_LOH=rbind(non_LOH,non_loh)
      }
    }
  } else {non_LOH=data.frame(start=chr_interval[1],end=chr_interval[2])} # in case no LOH is identified by IVD-PCF
  if (nrow(non_LOH)>0){
    for (j in 1:nrow(non_LOH)){
      if (non_LOH$start[j]<chr_loc[i,]$cen.left.base & non_LOH$end[j]>chr_loc[i,]$cen.right.base){
        start.pos=c(non_LOH$start[j],chr_loc[i,]$cen.right.base)
        end.pos=c(chr_loc[i,]$cen.left.base,non_LOH$end[j])
        non_LOH=non_LOH[-j,]
        non_LOH=rbind(non_LOH, data.frame(start=start.pos,end=end.pos))
      }
    }
    non_LOH$diff=non_LOH$end-non_LOH$start
  }
  
  non_LOH=non_LOH[order(non_LOH$start),] # the non_LOH should always be in order by position
  
  #STEP 2.1: identify LOH by inter-HET SNP regions # differentiating from HOM stretch in sample with logR < -0.8
  winsize=MIN_HET_DIST 
  ohet=GL_OHET[[i]] 
  nSNPs=as.numeric(nrow(GL_LogR))
  logr=GL_LogR[which(GL_LogR$Chromosome==i),]
  colnames(logr)[3]="LogR"
  logr$Position=as.numeric(logr$Position)
  if (!is.null(non_LOH)){ # if regions of non_LOH exist after IVD-PCF, run window-based search 
    pLOH_regions=data.frame()
    if (is.na(match(i,c(13,14,15,21,22)))){
      print(paste("START",i,"p ARM"))
      PARM=non_LOH[which(non_LOH$end<=chr_loc[i,]$cen.left.base),]
      if (nrow(PARM)>0){
        #if (nrow(PARM)==1 & non_LOH$start[1]==chr_interval[1] & non_LOH$end[1]==chr_interval[2]){
        parm=PARM
      } else if (nrow(PARM)==0 & sum(non_LOH$diff)!=0) {
        parm=data.frame(start=chr_interval[1],end=chr_loc[i,]$cen.left.base-CENTROMERE_DIST)
      } else {print("unknown issue")}
      
      if (parm[nrow(parm),1]<(parm[nrow(parm),2]-CENTROMERE_DIST)){ 
        parm[nrow(parm),2]=parm[nrow(parm),2]-CENTROMERE_DIST # exclude the last CENTROMERE_DIST segment next to the centromere (left side) - too noisy
      } else {parm=parm[-nrow(parm),]}
      parm$diff=parm$end-parm$start
      
      # search per non_LOH segment
      for (seg in 1:nrow(parm)){
        LoH=data.frame()
        #IVD-based breakpoints for small regions#
        seg_ivd=ohet[which(ohet$Position_dist>=MIN_HET_DIST & ohet$Position>=parm$start[seg] & ohet$Position<=parm$end[seg]),]
        #if (!is.null(nrow(seg_ivd))){
        if (nrow(seg_ivd)>0){
          win=nrow(seg_ivd)
          print(win)
          # win=floor(parm$diff[seg]/winsize)
          # print(win)
          #if (win>0){
          for (j in 1:win){
            loh=NULL
            start=seg_ivd$Position[j]
            end=start+seg_ivd$Position_dist[j]
            COV=logr[which(logr$Position>start & logr$Position<end),] # logR of homozygote SNPs within
            medcov=median(COV[,3])
            cov=mean(COV[,3])
            denSNP=nrow(COV)/(nSNPs/sum(chr_loc$length)*seg_ivd$Position_dist[j])
            #if (!is.na(cov) & cov < -0.8 & medcov < -0.8 & !is.null(denSNP) & denSNP>0.5){ # to use a minimum SNP density of 0.5 to get logR estimate
            if (!is.na(cov) & !is.null(denSNP) & denSNP>0.5){ # to use a minimum SNP density of 0.5 to get logR estimate AND not put the cov cut-off before applying PCF
              #loh=data.frame(start=start,end=end,LogR=cov,medianLogR=medcov,denSNP=denSNP)
              jpcf=pcf(COV,gamma=GAMMA_LOGR,verbose = F)
              jpcf=jpcf[which(jpcf$mean < -0.8),]
              if (nrow(jpcf)>0){
                loh=data.frame(start=jpcf$start.pos[1],end=jpcf$end.pos[nrow(jpcf)],LogR=mean(jpcf$mean),denSNP=denSNP)
                loh$N=nrow(logr[which(logr$Position>=loh$start & logr$Position<=loh$end),])
                if (loh$N<10){loh=NULL} # if LOH region is supported by less than 10 SNPs, then remove it
              }
            }
            if (!is.null(loh)){
              LoH=rbind(LoH,loh)
            }
            if (j %% 100 ==0){
              print(paste("interval=",j))
            }
          }
        } else {print(paste("no het SNPs in segment",seg))}
        # no. of LOH intervals
        print(paste("p-arm nrow(LOH) segment",seg,"=",nrow(LoH)))
        if (nrow(LoH)==0){
          print(paste("No LOH identified in p-arm segment",seg))
        } else{
          if (nrow(LoH)==1){
            LoH_regions=data.frame(chrom=i,arm="p",start.pos=LoH$start,end.pos=LoH$end)
          }
          if (nrow(LoH)>1){
            #combine smaller regions into larger regions of LOH
            LoH_regions=data.frame()
            start=LoH$start[1]
            for (j in 2:nrow(LoH)){
              print(j)
              if (LoH$start[j]==LoH$end[j-1]){
                end=LoH$end[j] # include the new row (i) in the merge
              }
              else {
                end=LoH$end[j-1] # stop merge at the previous row (i-1)
                LoH_regions=rbind(LoH_regions,data.frame(chrom=i,arm="p",start.pos=start,end.pos=end))
                start=LoH$start[j]
              }
            }
            # add final block if it ends at the end of the LoH dataframe
            if (end==LoH$end[nrow(LoH)]){
              LoH_regions=rbind(LoH_regions,data.frame(chrom=i,arm="p",start.pos=start,end.pos=end))
            }
            else if (start==LoH$start[nrow(LoH)] & end==LoH$end[nrow(LoH)-1]){
              LoH_regions=rbind(LoH_regions,data.frame(chrom=i,arm="p",start.pos=start,end.pos=LoH$end[nrow(LoH)]))
            }
          }
          pLOH_regions=rbind(pLOH_regions,LoH_regions)
        }
      }
      if (nrow(pLOH_regions)>0){
        #pARM BAF/LogR plot(s)
        pdf(paste0(GERMLINENAME,"_chr",i,"_",MIN_HET_DIST/1e3,"k_based_pLOH_events.pdf"))
        suppressWarnings(
          for (s in 1:nrow(pLOH_regions)){
            sBAF=ggplot(ohet,aes(Position,baf))+geom_jitter()+ylim(0,1)+
              geom_vline(xintercept = c(pLOH_regions$start.pos[s],pLOH_regions$end.pos[s]),col="red",linetype="longdash")+
              xlim(pLOH_regions$start.pos[s]-LENGTH_ADJACENT,pLOH_regions$end.pos[s]+LENGTH_ADJACENT)+
              ggtitle(paste("pARM LOH region",s))+labs(y="BAF")
            sLogR=ggplot(logr,aes(Position,LogR))+geom_jitter()+ylim(-5.2,1.2)+
              geom_vline(xintercept = c(pLOH_regions$start.pos[s],pLOH_regions$end.pos[s]),col="red",linetype="longdash")+
              xlim(pLOH_regions$start.pos[s]-LENGTH_ADJACENT,pLOH_regions$end.pos[s]+LENGTH_ADJACENT)
            grid.newpage()
            grid.draw(rbind(ggplotGrob(sBAF), ggplotGrob(sLogR), size = "last"))
            #print(plot_grid(sBAF,sLogR, ncol = 1, align = "v"))
          }
        )
        dev.off()
        #
        print("Candidate LOH regions plotted for pARM")
      }
    } else {print(paste("chr",i,"is acrocentric - no p arm analysis"))}
    # Q ARM RUN:
    print(paste("START",i,"q ARM"))
    qLOH_regions=data.frame()
    QARM=non_LOH[which(non_LOH$start>=chr_loc[i,]$cen.right.base),]
    if (nrow(QARM)>0){
      #if (nrow(PARM)==1 & non_LOH$start[1]==chr_interval[1] & non_LOH$end[1]==chr_interval[2]){
      qarm=QARM
    } else if (nrow(QARM)==0 & sum(non_LOH$diff)!=0) {
      qarm=data.frame(start=chr_loc[i,]$cen.right.base,end=chr_interval[2])
    } else {print("unknown issue")}
    qarm[1,1]=qarm[1,1]+CENTROMERE_DIST # to exclude the first CENTROMERE_DIST next to the centromere (right side) - noisy
    qarm$diff=qarm$end-qarm$start
    #
    # search per non_LOH segment
    for (seg in 1:nrow(qarm)){
      LoH=data.frame()
      #IVD-based breakpoints for small regions#
      seg_ivd=ohet[which(ohet$Position_dist>=MIN_HET_DIST & ohet$Position>=qarm$start[seg] & ohet$Position<=qarm$end[seg]),]
      #if (!is.null(nrow(seg_ivd))){
      if (nrow(seg_ivd)>0){
        win=nrow(seg_ivd)
        print(win)
        # win=floor(qarm$diff[seg]/winsize)
        # print(win)
        #if (win>0){
        for (j in 1:win){
          loh=NULL
          start=seg_ivd$Position[j]
          end=start+seg_ivd$Position_dist[j]
          COV=logr[which(logr$Position>start & logr$Position<end),] # logR of homozygote SNPs within
          cov=mean(COV[,3])
          medcov=median(COV[,3])
          denSNP=nrow(COV)/(nSNPs/sum(chr_loc$length)*seg_ivd$Position_dist[j])
          #if (!is.na(cov) & cov < -0.8 & medcov < -0.8 & !is.null(denSNP) & denSNP>0.5){ # to use a minimum SNP density of 0.5 to get logR estimate
           if (!is.na(cov) & !is.null(denSNP) & denSNP>0.5){ # to use a minimum SNP density of 0.5 to get logR estimate AND not put the cov cut-off before applying PCF
            #loh=data.frame(start=start,end=end,LogR=cov,medianLogR=medcov,denSNP=denSNP)
            jpcf=pcf(COV,gamma=GAMMA_LOGR,verbose = F)
            jpcf=jpcf[which(jpcf$mean < -0.8),]
            if (nrow(jpcf)>0){
              loh=data.frame(start=jpcf$start.pos[1],end=jpcf$end.pos[nrow(jpcf)],LogR=mean(jpcf$mean),denSNP=denSNP)
              loh$N=nrow(logr[which(logr$Position>=loh$start & logr$Position<=loh$end),])
              if (loh$N<10){loh=NULL} # if LOH region is supported by less than 10 SNPs, then remove it
            }
          }
          if (!is.null(loh)){
            LoH=rbind(LoH,loh)
          }
          if (j %% 100 ==0){
            print(paste("interval=",j))
          }
        }
      } else {print(paste("no het SNPs in segment",seg))}
      
      # no. of LoH intervals
      print(paste("q-arm nrow(LoH) segment",seg,"=",nrow(LoH)))
      if (nrow(LoH)==0){
        print(paste("No LOH identified in q-arm segment",seg))
      } else {
        if (nrow(LoH)==1){
          LoH_regions=data.frame(chrom=i,arm="q",start.pos=LoH$start,end.pos=LoH$end)
        }
        if (nrow(LoH)>1){
          LoH_regions=data.frame()
          #combine smaller regions into larger regions of LOH
          start=LoH$start[1]
          for (j in 2:nrow(LoH)){
            print(j)
            if (LoH$start[j]==LoH$end[j-1]){
              end=LoH$end[j] # include the new row (i) in the merge
            }
            else {
              end=LoH$end[j-1] # stop merge at the previous row (i-1)
              LoH_regions=rbind(LoH_regions,data.frame(chrom=i,arm="q",start.pos=start,end.pos=end))
              start=LoH$start[j]
            }
          }
          # add final block if it ends at the end of the LOH dataframe
          if (end==LoH$end[nrow(LoH)]){
            LoH_regions=rbind(LoH_regions,data.frame(chrom=i,arm="q",start.pos=start,end.pos=end))
          }
          else if (start==LoH$start[nrow(LoH)] & end==LoH$end[nrow(LoH)-1]){
            LoH_regions=rbind(LoH_regions,data.frame(chrom=i,arm="q",start.pos=start,end.pos=LoH$end[nrow(LoH)]))
          }
        }
        qLOH_regions=rbind(qLOH_regions,LoH_regions)
      }
    }
    if (nrow(qLOH_regions)>0){
      #qARM BAF/LogR plot(s)
      pdf(paste0(GERMLINENAME,"_chr",i,"_",MIN_HET_DIST/1e3,"k_based_qLOH_events.pdf"))
      suppressWarnings(
        for (s in 1:nrow(qLOH_regions)){
          sBAF=ggplot(ohet,aes(Position,baf))+geom_jitter()+ylim(0,1)+
            geom_vline(xintercept = c(qLOH_regions$start.pos[s],qLOH_regions$end.pos[s]),col="red",linetype="longdash")+
            xlim(qLOH_regions$start.pos[s]-LENGTH_ADJACENT,qLOH_regions$end.pos[s]+LENGTH_ADJACENT)+
            ggtitle(paste("qARM LOH region",s))
          sLogR=ggplot(logr,aes(Position,LogR))+geom_jitter()+ylim(-5.2,1.2)+
            geom_vline(xintercept = c(qLOH_regions$start.pos[s],qLOH_regions$end.pos[s]),col="red",linetype="longdash")+
            xlim(qLOH_regions$start.pos[s]-LENGTH_ADJACENT,qLOH_regions$end.pos[s]+LENGTH_ADJACENT)
          grid.newpage()
          grid.draw(rbind(ggplotGrob(sBAF), ggplotGrob(sLogR), size = "last"))
          #print(plot_grid(sBAF,sLogR, ncol = 1, align = "v"))
        }
      )
      dev.off()
      #
      print("Candidate LOH regions plotted for qARM")
    }
    #STEP 2.2: clean-up LOH[[i]] and merge LOH regions of both methods
    
    noLOH=NULL
    if (!is.null(nrow(LOH[[i]]))){
      for (j in 1:nrow(LOH[[i]])){
        LOH[[i]]$logR[j]=mean(logr[which(logr$Position>=LOH[[i]]$start.pos[j] & logr$Position<=LOH[[i]]$end.pos[j]),][,3])
        LOH[[i]]$nSNP[j]=nrow(logr[which(logr$Position>=LOH[[i]]$start.pos[j] & logr$Position<=LOH[[i]]$end.pos[j]),])
        LOH[[i]]$denSNP[j]=LOH[[i]]$nSNP[j]/((LOH[[i]]$end.pos[j]-LOH[[i]]$start.pos[j])*nSNPs/sum(chr_loc$length))
        if (LOH[[i]]$logR[j]>-0.8 | LOH[[i]]$denSNP[j]<0.5){
          noLOH=append(noLOH,j)
          print(j)
        }
      }
      LOH[[i]]=LOH[[i]][-noLOH,]
      LOH[[i]]=ifelse(nrow(LOH[[i]])==0,0,LOH[[i]])
    }
    
    LOH_regions=data.frame()
    if (nrow(pLOH_regions)>0){
      LOH_regions=rbind(LOH_regions,pLOH_regions)
    } else {print("no window-based LOH regions identified in p arm of non_LOH of IVD-PCF")}
    if (nrow(qLOH_regions)>0){
      LOH_regions=rbind(LOH_regions,qLOH_regions)
    } else {print("no window-based LOH regions identified in q arm of non_LOH of IVD-PCF")}
    if (nrow(LOH_regions)>0){
      if (!is.null(nrow(LOH[[i]]))){
        LOH[[i]]=rbind(LOH[[i]][,c("chrom","arm","start.pos","end.pos")],LOH_regions)
        LOH[[i]]=LOH[[i]][order(LOH[[i]]$start.pos),]
      } else {
        LOH[[i]]=LOH_regions
      }
    }
    
    
    #combine adjacent regions into larger regions of LOH
    if (!is.null(nrow(LOH[[i]]))){
      LOH[[i]]=LOH[[i]][!duplicated(LOH[[i]]),]
      LOHall=data.frame()
      ChrArms=unique(LOH[[i]]$arm)
      for (arm in ChrArms){
        LOHarm=LOH[[i]][LOH[[i]]$arm==arm,]
        if (nrow(LOHarm)>1){
          start=LOHarm$start.pos[1]
          for (j in 2:nrow(LOHarm)){
            print(j)
            if (LOHarm$start.pos[j]==LOHarm$end.pos[j-1]){
              end=LOHarm$end.pos[j] # include the new row (i) in the merge
            } else {
              if (LOHarm$start.pos[j]>LOHarm$end.pos[j-1]){
                end=LOHarm$end.pos[j-1] # stop merge at the previous row (i-1)
                LOHall=rbind(LOHall,data.frame(chrom=i,arm=arm,start.pos=start,end.pos=end))
                start=LOHarm$start.pos[j]
              } else if (LOHarm$start.pos[j]<LOHarm$end.pos[j-1]){
                end=max(LOHarm$end.pos[j-1],LOHarm$end.pos[j])
                start=min(start,LOHarm$start.pos[j])
                LOHall=rbind(LOHall,data.frame(chrom=i,arm=arm,start.pos=start,end.pos=end))
              }
            }
          }
          # add final block if it ends at the end of the LoH dataframe
          if (end==LOHarm$end.pos[nrow(LOHarm)]){
            LOHall=rbind(LOHall,data.frame(chrom=i,arm=arm,start.pos=start,end.pos=end))
          } else if (start==LOHarm$start.pos[nrow(LOHarm)] & end==LOHarm$end.pos[nrow(LOHarm)-1]){
            LOHall=rbind(LOHall,data.frame(chrom=i,arm=arm,start.pos=start,end.pos=LOHarm$end.pos[nrow(LOHarm)]))
          }
        } else {LOHall=rbind(LOHall,LOHarm[,c("chrom","arm","start.pos","end.pos")])}
      }
    } else {LOHall=LOH[[i]]}
    print("LOHall")
    print(LOHall)
  } else { # no non_LOH region was found - all chromosome is called as LOH (highly unlikely at germline level)
    LOHall=LOH[[i]][,c("chrom","arm","start.pos","end.pos")]
    print("LOHall")
    print(LOHall)
  }
  
  if (!is.null(nrow(LOHall))){
    LOHall=LOHall[!duplicated(LOHall),]
    LOHall$diff=LOHall$end.pos-LOHall$start.pos
  } else {print(paste("no LOH (IVD and/or window-based) was identified for chr",i))}
  if (exists("non_loh")){
    rm(non_loh)}
  if (exists("non_LOH")){
    rm(non_LOH)
  }
  
  #STEP 3####################################################################################################################################################
  # RECONSTRUCT alleleCounter files for the pseudo-NORMAL sample
  # use loop to find intervening blocks with no LOH - while taking account of the centromere - RUN2#
  if (!is.null(nrow(LOHall))){
    names(ac)=c("chr","position",1:4,"depth")
    chr_interval=c(ac$position[1],ac$position[nrow(ac)])
    non_LOH=data.frame()####################################### get all non_LOH regions#
    for (j in 1:(nrow(LOHall)+1)){
      if (j == 1 & chr_interval[1]==LOHall$start.pos[j]){
        print("LOH from start of chromosome")
      } else if (j == 1 & chr_interval[1]<LOHall$start.pos[j]){
        non_loh=data.frame(start=chr_interval[1],end=LOHall$start.pos[j]-1)
        print("ONE")
      } else if (j>1 & j <= nrow(LOHall) & LOHall$arm[j]==LOHall$arm[j-1]){
        non_loh=data.frame(start=LOHall$end.pos[j-1]+1,end=LOHall$start.pos[j]-1)
        print("TWO")
      } else if (j>1 & j <= nrow(LOHall) & LOHall$arm[j]!=LOHall$arm[j-1]){
        non_loh=data.frame(start=c(min(LOHall$end.pos[j-1]+1,chr_loc[i,]$cen.left.base),chr_loc[i,]$cen.right.base),end=c(chr_loc[i,]$cen.left.base,LOHall$start.pos[j]-1))
        print("THREE")
      } else{
        if ((LOHall$end.pos[j-1]+1)<chr_interval[2]){ # avoids going over the chromosome interval
          non_loh=data.frame(start=LOHall$end.pos[j-1]+1,end=chr_interval[2])
        } else{
          print("reached end of chromosome")
          rm(non_loh)
        }
      }
      print(j)
      if (exists("non_loh")){
        non_LOH=rbind(non_LOH,non_loh)
      }
    } 
    # the non-LOH region length from PCF is:
    if (!is.null(nrow(non_LOH))){
      non_LOH$length=non_LOH$end-non_LOH$start
      non_LOH=non_LOH[non_LOH$length>=0,] # picks all non_LOH segments even if 1bp in length
      non_LOH_length=sum(non_LOH$length) # total length of non-LOH regions in chr i
      print(paste("Total length of non LOH regions =",non_LOH_length))
      # average Het SNP interval:
      if (non_LOH_length>1e6){ # run this only if combined non-LOH regions are at least 1Mb long
        SNP_interval=non_LOH_length/nrow(GL_OHET[[i]]) # estimate of genomic space between any two Het SNPs
      } else {SNP_interval = 2000} # replace with 5000 to increase run speed!?
      # no. of SNPs to be Hets in the LOH region (COMBINED FOR THE WHOLE CHROMOSOME):
      LOH_hetSNP_number=floor(sum(LOHall$diff)/SNP_interval)
      print(paste("No. of Het SNPs to be added to LOH regions:",LOH_hetSNP_number))
    }
    # reconstruct allele counts for the LOH region based on actual depth for all to be perfect heterozygotes - allele counts remain as integers
    lohs=data.frame() ####################################### get all non_LOH regions####
    for (j in 1:nrow(LOHall)){
      loh=ac[which(ac$position>=LOHall$start.pos[j] & ac$position<=LOHall$end.pos[j]),]
      m=merge(loh,al,"position")
      if (nrow(m)==nrow(loh)){
        print("merge OK")
      } else {print("ERROR - merge not OK")}
      # RE-reconstruct allele counts for LOH region
      hetSNP_number=max(LOHall$diff[j]/SNP_interval,10) #at least ten SNPs (if available in region) should be spiked in to be heterozygotes for PCF in Battenberg to pick it up
      if (nrow(m)>=hetSNP_number){
        print("more rows in LOH region than Het SNP number")
        spike=c(1,head(which(1:nrow(m) %% floor(nrow(m)/(hetSNP_number-1))==0),-1),nrow(m)) # to make the exact breakpoints are seen by Battenberg - making 1st and last SNP in region heterozygote
        for (k in spike){
          #for (k in 1:nrow(m)){
          #if (k %% floor(nrow(m)/hetSNP_number)==0){
          m$depth[k]=max(m$depth[k],10)
          m[cbind(k,2+m$a0[k])]=ifelse(m$depth[k]%%2==0,m$depth[k]/2,ceiling(m$depth[k]/2))
          m[cbind(k,2+m$a1[k])]=ifelse(m$depth[k]%%2==0,m$depth[k]/2,floor(m$depth[k]/2))
          print(k)
          #}
        }
      } else {
        print("less rows in LOH region than Het SNP number - turning all into Heterozygotes") # technically shouldn't happen 
        for (k in 1:nrow(m)){
          m$depth[k]=max(m$depth[k],10)
          m[cbind(k,2+m$a0[k])]=ifelse(m$depth[k]%%2==0,m$depth[k]/2,ceiling(m$depth[k]/2))
          m[cbind(k,2+m$a1[k])]=ifelse(m$depth[k]%%2==0,m$depth[k]/2,floor(m$depth[k]/2))
          print(k)
        }
      }
      print(paste("LOH region segment",j))
      lohs=rbind(lohs,m)
    }
    
    lohs=lohs[,c("chr","position",1:4,"depth")]
    ####
    # combine alleleCounts for LOHS and non_LOH regions####
    non_lohs=data.frame()
    for (j in 1:nrow(non_LOH)){
      non_loh=ac[which(ac$position>=non_LOH$start[j] & ac$position<=non_LOH$end[j]),]
      non_lohs=rbind(non_lohs,non_loh)
      print(paste("non_LOH segment",j,"added"))
    }
    #write out as alleleCounts file - "normal" ID #####################################
    if (nrow(non_lohs)+nrow(lohs)==nrow(ac)){
      ac_out=rbind(non_lohs,lohs)
      ac_out=ac_out[order(ac_out$position),]
      write.table(ac_out,paste0(NORMALNAME,"_alleleFrequencies_chr",i,".txt"),col.names=F,row.names=F,quote=F,sep="\t")
      print(paste("reconstruction OK - new alleleCounts file generated for chr",i))
    } else {
      centro_ac=ac[which(ac$position>chr_loc$cen.left.base[i] & ac$position<chr_loc$cen.right.base[i]),]
      ac_out=rbind(non_lohs,lohs,centro_ac)
      ac_out=ac_out[order(ac_out$position),]
      ac_out=ac_out[!duplicated(ac_out$position),]
      if (nrow(ac_out)==nrow(ac)){
        print("reconstruction OK but SNPs found in the centromeric region - adding them back for consistency with original ac files")
        write.table(ac_out,paste0(NORMALNAME,"_alleleFrequencies_chr",i,".txt"),col.names=F,row.names=F,quote=F,sep="\t")
        } else {
      print("ERROR - missing SNPs - LOH and non-LOH regions not generated correctly; no AC file generated")}
      }
  } else {
    ac_out=ac
    write.table(ac_out,paste0(NORMALNAME,"_alleleFrequencies_chr",i,".txt"),col.names=F,row.names=F,quote=F,sep="\t")
    print(paste("No changes made to the alleleCounter file - no LOH in chr",i))
  }
  print(paste("STEP 2&3 - chr",i,"completed"))
}

#' Prepare data for impute
#'
#' @param chrom The chromosome for which impute input should be generated.
#' @param germline.allele.counts.file Output from the allele counter on the matched germline for this chromosome.
#' @param normal.allele.counts.file Output from the allele counter on the matched normal for this chromosome.
#' @param output.file File where the impute input for this chromosome will be written.
#' @param imputeinfofile Info file with impute reference information.
#' @param is.male Boolean denoting whether this sample is male (TRUE), or female (FALSE).
#' @param problemLociFile A file containing genomic locations that must be discarded (optional).
#' @param useLociFile A file containing genomic locations that must be included (optional).
#' @param heterozygousFilter The cutoff where a SNP will be considered as heterozygous (default 0.01).
#' @author dw9, sd11, Naser Ansari-Pour (BDI, Oxford)
#' @export
generate.impute.input.wgs.germline = function(chrom, germline.allele.counts.file, normal.allele.counts.file, output.file, imputeinfofile, is.male, problemLociFile=NA, useLociFile=NA, heterozygousFilter=0.1) {
  
  # Read in the 1000 genomes reference file paths for the specified chrom
  impute.info = parse.imputeinfofile(imputeinfofile, is.male, chrom=chrom)
  chr_names = unique(impute.info$chrom)
  chrom_name = parse.imputeinfofile(imputeinfofile, is.male)$chrom[chrom]
  
  #print(paste("GenerateImputeInput is.male? ", is.male,sep=""))
  #print(paste("GenerateImputeInput #impute files? ", nrow(impute.info),sep=""))
  
  # Read in the known SNP locations from the 1000 genomes reference files
  known_SNPs = read.table(impute.info$impute_legend[1], sep=" ", header=T, stringsAsFactors=F)
  if(nrow(impute.info)>1){
    for(r in 2:nrow(impute.info)){
      known_SNPs = rbind(known_SNPs, read.table(impute.info$impute_legend[r], sep=" ", header=T, stringsAsFactors=F))
    }
  }
  
  # filter out bad SNPs (streaks in BAF)
  if((problemLociFile != "NA") & (!is.na(problemLociFile))) {
    problemSNPs = read.table(problemLociFile, header=T, sep="\t", stringsAsFactors=F)
    problemSNPs = problemSNPs$Pos[problemSNPs$Chr==chrom_name]
    badIndices = match(known_SNPs$position, problemSNPs)
    known_SNPs = known_SNPs[is.na(badIndices),]
    rm(problemSNPs, badIndices)
  }
  
  # filter 'good' SNPs (e.g. SNP6 positions)
  if((useLociFile != "NA") & (!is.na(useLociFile))) {
    goodSNPs = read.table(useLociFile, header=T, sep="\t", stringsAsFactors=F)
    goodSNPs = goodSNPs$pos[goodSNPs$chr==chrom_name]  
    len = length(goodSNPs)
    goodIndices = match(known_SNPs$position, goodSNPs)
    known_SNPs = known_SNPs[!is.na(goodIndices),]
    rm(goodSNPs, goodIndices)
  }
  
  # Read in the allele counts and see which known SNPs are covered
  snp_data = read.table(germline.allele.counts.file, comment.char="#", sep="\t", header=F, stringsAsFactors=F)
  normal_snp_data = read.table(normal.allele.counts.file, comment.char="#", sep="\t", header=F, stringsAsFactors=F)
  snp_data = cbind(snp_data, normal_snp_data)
  indices = match(known_SNPs$position, snp_data[,2])
  found_snp_data = snp_data[indices[!is.na(indices)],]
  rm(snp_data)
  
  # Obtain BAF for this chromosome (note: this is quicker than reading in the whole genome BAF file generated in the earlier step)
  nucleotides = c("A","C","G","T")
  ref_indices = match(known_SNPs[!is.na(indices),3], nucleotides)+ncol(normal_snp_data)+2
  alt_indices = match(known_SNPs[!is.na(indices),4], nucleotides)+ncol(normal_snp_data)+2
  BAFs = as.numeric(found_snp_data[cbind(1:nrow(found_snp_data),alt_indices)])/(as.numeric(found_snp_data[cbind(1:nrow(found_snp_data),alt_indices)])+as.numeric(found_snp_data[cbind(1:nrow(found_snp_data),ref_indices)]))
  BAFs[is.nan(BAFs)] = 0
  rm(nucleotides, ref_indices, alt_indices, found_snp_data, normal_snp_data)
  
  # Set the minimum level to use for obtaining genotypes
  minBaf = min(heterozygousFilter, 1.0-heterozygousFilter)
  maxBaf = max(heterozygousFilter, 1.0-heterozygousFilter)
  
  # Obtain genotypes that impute2 is able to understand
  genotypes = array(0,c(sum(!is.na(indices)),3))
  genotypes[BAFs<=minBaf,1] = 1
  genotypes[BAFs>minBaf & BAFs<maxBaf,2] = 1
  genotypes[BAFs>=maxBaf,3] = 1	
  
  # Create the output
  snp.names = paste("snp",1:sum(!is.na(indices)), sep="")
  out.data = cbind(snp.names, known_SNPs[!is.na(indices),1:4], genotypes)
  
  write.table(out.data, file=output.file, row.names=F, col.names=F, quote=F)
  if(is.na(as.numeric(chrom_name))) {
    sample.g.file = paste(dirname(output.file), "/sample_g.txt", sep="")
    #not sure this is necessary, because only the PAR regions are used for males
    #if(is.male){
    #	sample_g_data=data.frame(ID_1=c(0,"INDIVI1"),ID_2=c(0,"INDIVI1"),missing=c(0,0),sex=c("D",1))
    #}else{
    sample_g_data = data.frame(ID_1=c(0,"INDIVI1"), ID_2=c(0,"INDIVI1"), missing=c(0,0), sex=c("D",2))
    #}
    write.table(sample_g_data, file=sample.g.file, row.names=F, col.names=T, quote=F)
  }
}

#' Function to correct LogR for waivyness that correlates with GC content
#' @param germline_LogR_file String pointing to the germline LogR output
#' @param outfile String pointing to where the GC corrected LogR should be written
#' @param correlations_outfile File where correlations are to be saved
#' @param gc_content_file_prefix String pointing to where GC windows for this reference genome can be 
#' found. These files should be split per chromosome and this prefix must contain the full path until
#' chr in its name. The .txt extension is automatically added.
#' @param replic_timing_file_prefix Like the gc_content_file_prefix, containing replication timing info (supply NULL if no replication timing correction is to be applied)
#' @param chrom_names A vector containing chromosome names to be considered
#' @param recalc_corr_afterwards Set to TRUE to recalculate correlations after correction
#' @author jonas demeulemeester, sd11, Naser Ansari-Pour (BDI, Oxford)
#' @export
gc.correct.wgs.germline = function(germline_LogR_file, outfile, correlations_outfile, gc_content_file_prefix, replic_timing_file_prefix, chrom_names, recalc_corr_afterwards=F) {
  
  if (is.null(gc_content_file_prefix)) {
    stop("GC content reference files must be supplied to WGS GC content correction")
  }
  
  Germline_LogR = read_logr(germline_LogR_file)
  
  print("Processing GC content data")
  chrom_idx = 1:length(chrom_names)
  gc_files = paste0(gc_content_file_prefix, chrom_idx, ".txt.gz")
  GC_data = do.call(rbind, lapply(gc_files, read_gccontent))
  colnames(GC_data) = c("chr", "Position", paste0(c(25,50,100,200,500), "bp"),
                        paste0(c(1,2,5,10,20,50,100), "kb"))#,200,500), "kb"),
  # paste0(c(1,2,5,10), "Mb"))
  
  if (!is.null(replic_timing_file_prefix)) {
    print("Processing replciation timing data")
    replic_files = paste0(replic_timing_file_prefix, chrom_idx, ".txt.gz")
    replic_data = do.call(rbind, lapply(replic_files, read_replication))
  }
  
  # omit non-matching loci, replication data generated at exactly same GC loci
  locimatches = match(x = paste0(Germline_LogR$Chromosome, "_", Germline_LogR$Position),
                      table = paste0(GC_data$chr, "_", GC_data$Position))
  Germline_LogR = Germline_LogR[which(!is.na(locimatches)), ]
  GC_data = GC_data[na.omit(locimatches), ]
  if (!is.null(replic_timing_file_prefix)) {
    replic_data = replic_data[na.omit(locimatches), ]
  }
  rm(locimatches)
  
  corr = abs(cor(GC_data[, 3:ncol(GC_data)], Germline_LogR[,3], use="complete.obs")[,1])
  if (!is.null(replic_timing_file_prefix)) {
    corr_rep = abs(cor(replic_data[, 3:ncol(replic_data)], Germline_LogR[,3], use="complete.obs")[,1])
  }
  
  index_1kb = which(names(corr)=="1kb")
  maxGCcol_insert = names(which.max(corr[1:index_1kb]))
  index_100kb = which(names(corr)=="100kb")
  # start large window sizes at 5kb rather than 2kb to avoid overly correlated expl variables
  maxGCcol_amplic = names(which.max(corr[(index_1kb+2):index_100kb]))
  if (!is.null(replic_timing_file_prefix)) {
    maxreplic = names(which.max(corr_rep))
  }
  
  if (!is.null(replic_timing_file_prefix)) {
    cat("Replication timing correlation: ",paste(names(corr_rep),format(corr_rep,digits=2), ";"),"\n") 
    cat("Replication dataset: " ,maxreplic,"\n")
  }
  cat("GC correlation: ",paste(names(corr),format(corr,digits=2), ";"),"\n")   
  cat("Short window size: ",maxGCcol_insert,"\n")
  cat("Long window size: ",maxGCcol_amplic,"\n")
  
  if (!is.null(replic_timing_file_prefix)) {
    # Multiple regression - with replication timing
    corrdata = data.frame(logr = Germline_LogR[,3, drop = T],
                          GC_insert = GC_data[,maxGCcol_insert, drop = T],
                          GC_amplic = GC_data[,maxGCcol_amplic, drop = T],
                          replic = replic_data[, maxreplic, drop = T])
    colnames(corrdata) = c("logr", "GC_insert", "GC_amplic", "replic")
    if (!recalc_corr_afterwards)
      rm(GC_data, replic_data)
    
    model = lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = T) + splines::ns(x = GC_amplic, df = 5, intercept = T) + splines::ns(x = replic, df = 5, intercept = T), y=F, model = F, data = corrdata, na.action="na.exclude")
    
    corr = data.frame(windowsize=c(names(corr), names(corr_rep)), correlation=c(corr, corr_rep))
    write.table(corr, file=gsub(".txt", "_beforeCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)
    
  } else {
    # Multiple regression  - without replication timing
    corrdata = data.frame(logr = Germline_LogR[,3, drop = T],
                          GC_insert = GC_data[,maxGCcol_insert, drop = T],
                          GC_amplic = GC_data[,maxGCcol_amplic, drop = T])
    colnames(corrdata) = c("logr", "GC_insert", "GC_amplic")
    if (!recalc_corr_afterwards)
      rm(GC_data)
    
    model = lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = T) + splines::ns(x = GC_amplic, df = 5, intercept = T), y=F, model = F, data = corrdata, na.action="na.exclude")
    
    corr = data.frame(windowsize=names(corr), correlation=corr)
    write.table(corr, file=gsub(".txt", "_beforeCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)
  }
  
  Germline_LogR[,3] = residuals(model)
  rm(model, corrdata)
  
  readr::write_tsv(x=Germline_LogR[which(!is.na(Germline_LogR[,3])), ], file=outfile)
  
  if (recalc_corr_afterwards) {
    # Recalculate the correlations to see how much there is left
    corr = abs(cor(GC_data[, 3:ncol(GC_data)], Germline_LogR[,3], use="complete.obs")[,1])
    if (!is.null(replic_timing_file_prefix)) {
      corr_rep = abs(cor(replic_data[, 3:ncol(replic_data)], Germline_LogR[,3], use="complete.obs")[,1])
      cat("Replication timing correlation post correction: ",paste(names(corr_rep),format(corr_rep,digits=2), ";"),"\n") 
    }
    cat("GC correlation post correction: ",paste(names(corr),format(corr,digits=2), ";"),"\n")   
    
    if (!is.null(replic_timing_file_prefix)) {
      corr = data.frame(windowsize=c(names(corr), names(corr_rep)), correlation=c(corr, corr_rep))
      write.table(corr, file=gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)
    } else {
      corr = data.frame(windowsize=c(names(corr)), correlation=corr)
      write.table(corr, file=gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)
    }
  } else {
    corr$correlation = NA
    write.table(corr, file=gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)
  }
}


#' Prepare WGS data of germline for haplotype construction
#' 
#' This function performs part of the Battenberg WGS pipeline: Counting alleles, generating BAF and logR, 
#' reconstructing normal-pair allele counts for the germline and performing GC content correction.
#' 
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param chrom_coord Full path to the file with chromosome coordinates including start, end and left/right centromere positions
#' @param germlinebam Full path to the germline BAM file 
#' @param germlinename Identifier to be used for germline output files (i.e. the germline BAM file name without the '.bam' extension).
#' @param g1000lociprefix Prefix path to the 1000 Genomes loci reference files
#' @param g1000allelesprefix Prefix path to the 1000 Genomes SNP allele reference files
#' @param gamma_ivd The PCF gamma value for segmentation of 1000G hetSNP IVD values (Default 1e5).
#' @param kmin_ivd The min number of SNPs to support a segment in PCF of 1000G hetSNP IVD values (Default 50)
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
prepare_wgs_germline = function(chrom_names, chrom_coord, germlinebam, germlinename, g1000lociprefix, g1000allelesprefix, gamma_ivd=1e5, kmin_ivd=50, centromere_noise_seg_size=1e6,
                                centromere_dist=5e5, min_het_dist=2e3, gamma_logr=100, length_adjacent=5e4, gccorrectprefix,repliccorrectprefix, min_base_qual, min_map_qual, 
                                 allelecounter_exe, min_normal_depth, skip_allele_counting) {
  
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  
  if (!skip_allele_counting) {
    # Obtain allele counts for 1000 Genomes locations for the germline
    foreach::foreach(i=1:length(chrom_names)) %dopar% {
      getAlleleCounts(bam.file=germlinebam,
                      output.file=paste(germlinename,"_alleleFrequencies_chr", i, ".txt", sep=""),
                      g1000.loci=paste(g1000lociprefix, i, ".txt", sep=""),
                      min.base.qual=min_base_qual,
                      min.map.qual=min_map_qual,
                      allelecounter.exe=allelecounter_exe)
    }
  }
  
  # Standardise Chr notation (removes 'chr' string if present; essential for cell_line_baf_logR)

  standardiseChrNotation_germline(GERMLINENAME=germlinename)
  
  # Obtain BAF and LogR from the raw allele counts of the germline
  germline_baf_logR(GERMLINENAME=germlinename,
                     g1000alleles.prefix=g1000allelesprefix,
                     chrom_names=chrom_names
  )
  
  # Reconstruct normal-pair allele count files for the germline
  
  foreach::foreach(i=1:length(chrom_names),.export=c("germline_reconstruct_normal","GL_OHET","GL_AL","GL_AC","GL_LogR"),.packages=c("copynumber","ggplot2","grid")) %dopar% {
    
     germline_reconstruct_normal(GERMLINENAME=germlinename,
                                 NORMALNAME=paste0(germlinename,"_normal"),
                                 chrom_coord=chrom_coord,
                                 chrom=i,
                                 GL_OHET=GL_OHET,
                                 GL_AL=GL_AL,
                                 GL_AC=GL_AC,
                                 GL_LogR=GL_LogR,
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
  
  # Perform GC correction
  gc.correct.wgs.germline(germline_LogR_file=paste(germlinename,"_mutantLogR.tab", sep=""),
                 outfile=paste(germlinename,"_mutantLogR_gcCorrected.tab", sep=""),
                 correlations_outfile=paste(germlinename, "_GCwindowCorrelations.txt", sep=""),
                 gc_content_file_prefix=gccorrectprefix,
                 replic_timing_file_prefix=repliccorrectprefix,
                 chrom_names=chrom_names)
}
