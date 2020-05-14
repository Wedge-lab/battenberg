cell_line_reconstruct_normal <-function(TUMOURNAME,NORMALNAME,GAMMA,KMIN,CENTROMERE_DIST,MIN_HET_DIST,GAMMA_LOGR,LENGTH_ADJACENT){
  # IDENTIFY REGIONS OF LOH ####
  ##########################################################chrom=1:22
  colClasses=c(chr="numeric",start="numeric",cen.left.base="numeric",cen.right.base="numeric",end="numeric")
  chr_loc=read.table(paste0(Ref_files_dir,"gcCorrect_chromosome_coordinates_hg19.txt"),colClasses = colClasses,header=T,stringsAsFactors = F)
  chr_loc$length=(chr_loc$cen.left.base-chr_loc$start)+(chr_loc$end-chr_loc$cen.right.base)
  #STEP 2.0: identify LOH by IVD-PCF
  LOH=list()
  PCF_folder = "PCF_plots"
  if(!file.exists(PCF_folder)){
    dir.create(PCF_folder)
  }
  
  
  print(paste("chrom=",i))
  pcf_input=data.frame(chr=i,position=CL_OHET[[i]]$Position,IVD=(CL_OHET[[i]]$Position_dist_percent))
  pcf_input=pcf_input[which(pcf_input$position>=chr_loc[i,"start"] & pcf_input$position<=chr_loc[i,"end"]),] # use only regions covered with gcCorrect LogR range 
  PCF=pcf(pcf_input,gamma=GAMMA,kmin = KMIN)
  pdf(paste0(PCF_folder,"/",TUMOURNAME,"_chr",i,"_PCF_plot.pdf"))
  plotChrom(pcf_input,PCF)
  dev.off()
  PCF$diff=PCF$end.pos-PCF$start.pos
  
  # Decide if there is any LOH based on PCF and chr_snp_density
  chr_snp_density=nrow(pcf_input)/(pcf_input$position[nrow(pcf_input)]-pcf_input$position[1]) # density of HET SNPs across the region covered by HET SNPs
  ####CALCULATE min_normal_snp_density####
  # setwd("~/Documents/BDI/Oxford-Celgene_MM/")
  # l=list.files(pattern="hetero")
  # snp_per_chr=data.frame()
  # for (i in 1:length(l)){
  #   h=read.table(l[i],header=T,stringsAsFactors = F)
  #   hout=data.frame(chr=h$Chromosome[1],dist=h$Position[nrow(h)]-h$Position[1],no_snp=nrow(h))
  #   snp_per_chr=rbind(snp_per_chr,hout)
  #   print(h$Chromosome[1])
  # }
  # snp_per_chr$density=snp_per_chr$no_snp/snp_per_chr$dist
  # minimum normal density for SNPs (in bps) is 3 x 10^-4 with median of 7 x 10^-4
  ####
  min_normal_snp_density=0.0001
  loh_regions=PCF[which(round(PCF$mean,3)>0.001),] # LOH regions
  loh_regions=loh_regions[which(loh_regions$n.probes>1),] # only keep segments with minimum of 2 probes (SNPs) in PCF jump
  if (mean(pcf_input$IVD)>0.01 & chr_snp_density<min_normal_snp_density){ #can change chr_snp_density from 0.00005 to 0.0001 as conservative measure - done
    #mean(pcf_input$IVD) or mean(PCF$mean) indicates presence of jumps in IVD
    loh_regions=loh_regions # LOH regions
    print(paste("full-length chromosomal loss at chr",i))
  } else if (sum(loh_regions$diff)>=((pcf_input$position[nrow(pcf_input)]-pcf_input$position[1]))*0.9 & chr_snp_density>min_normal_snp_density){
    # do PCF regions cover >=90% of the chromosome & is the chromosome snp density above the minimum (see calculation above commented)
    loh_regions=0 # LOH regions
    print(paste("no PCF jumps at chr",i))
  } else {
    loh_regions=loh_regions # LOH regions
    print(paste("likely partial LOH(s) at chr",i))
  }
  
  # 02-01-2020 UPDATE: new loop to turn empty dataframe to 0 for loh_regions 
  suppressWarnings(
    if (loh_regions[1]!=0){
      if (nrow(loh_regions)==0){
        loh_regions=0
      } else {print("dataframe non-empty")}
    } else {print("no LOH at all")})
  
  #filter regions for those next to the centromere and 'short' (<200k)
  noise=NULL
  if (!is.null(nrow(loh_regions))){
    for (j in 1:nrow(loh_regions)){
      if (loh_regions$arm[j]=="p"){
        if (loh_regions$end.pos[j]-chr_loc$cen.left.base[i]<1e5 & loh_regions$diff[j]<1e6){ #FOR EXCLUSION: max distance to centromere = 100kb , max length of short LOH region = 1Mb
          noise=append(noise,j)
        }
      }
      if (loh_regions$arm[j]=="q"){
        if (loh_regions$start.pos[j]-chr_loc$cen.right.base[i]<1e5 & loh_regions$diff[j]<1e6){ #FOR EXCLUSION: max distance to centromere = 100kb , max length of short LOH region = 1Mb
          noise=append(noise,j)
        }
      }
    }
  } else {print("no 'centromere noise' calculation")}
  if (!is.null(noise)){
    LOH_regions=loh_regions[-noise,]
  } else {LOH_regions=loh_regions}
  ####
  #remove LOH regions in the p arm of acrocentric chromosomes 13,14,15,21 and 22
  if (!is.na(match(i,c(13:15,21:22))) & !is.null(nrow(LOH_regions))){
    LOH_regions=LOH_regions[which(LOH_regions$arm!="p"),]
  }
  ####
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
  ac=CL_AC[[i]]
  al=CL_AL[[i]]
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
  
  #STEP 2.1: identify LOH by windows (size=1e5)
  winsize=MIN_HET_DIST # optimum value is 1e5 in differentiating from HOM stretch in sample
  ohet=CL_OHET[[i]] 
  nSNPs=as.numeric(nrow(CL_LogR))
  logr=CL_LogR[which(CL_LogR$Chromosome==i),]
  colnames(logr)[3]="LogR"
  logr$Position=as.numeric(logr$Position)
  if (!is.null(non_LOH)){
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
        parm[nrow(parm),2]=parm[nrow(parm),2]-CENTROMERE_DIST # 30-03-2020 UPDATE: to exclude the last 500k next to the centromere (left side) - too noisy
      } else {parm=parm[-nrow(parm),]}
      #
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
            if (!is.na(cov) & cov < -0.8 & medcov < -0.8 & !is.null(denSNP) & denSNP>0.5){ # to use a minimum SNP density of 0.5 to get logR estimate #CLcode
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
        pdf(paste0(TUMOURNAME,"_chr",i,"_",MIN_HET_DIST/1e3,"k_based_pLOH_events.pdf"))
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
    qarm[1,1]=qarm[1,1]+CENTROMERE_DIST # 30-03-2020 UPDATE: to exclude the first 500k next to the centromere (right side) - noisy
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
          if (!is.na(cov) & cov < -0.8 & medcov < -0.8 & !is.null(denSNP) & denSNP>0.5){ # to use a minimum SNP density of 0.5 to get logR estimate #CLcode
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
      pdf(paste0(TUMOURNAME,"_chr",i,"_",MIN_HET_DIST/1e3,"k_based_qLOH_events.pdf"))
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
    #STEP 2.2: merge LOH regions of both methods
    LOH_regions=data.frame()
    if (nrow(pLOH_regions)>0){
      print(pLOH_regions)
      LOH_regions=rbind(LOH_regions,pLOH_regions)
    } else {print("no window-based LOH regions identified in p arm of non_LOH of IVD-PCF")}
    if (nrow(qLOH_regions)>0){
      print(qLOH_regions)
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
  } else { # no non_LOH region was found - all chromosome is called as LOH
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
  # use loop to find intervening blocks with no LOH - while taking account of the centromere - RUN2####
  if (!is.null(nrow(LOHall))){
    names(ac)=c("chr","position",1:4,"depth")
    chr_interval=c(ac$position[1],ac$position[nrow(ac)])
    non_LOH=data.frame()####################################### get all non_LOH regions#### 
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
      non_LOH=non_LOH[non_LOH$length>0,]
      non_LOH_length=sum(non_LOH$length) # total length of non-LOH regions in chr i
      print(paste("Total length of non LOH regions =",non_LOH_length))
      # average Het SNP interval:
      if (non_LOH_length>1e6){ # run this only if combined non-LOH regions are at least 1Mb long
        SNP_interval=non_LOH_length/nrow(CL_OHET[[i]]) # estimate of genomic space between any two Het SNPs
      } else {SNP_interval = 2000} # replace with 5000 to increase run speed!?
      # no. of SNPs to be Hets in the LOH region (COMBINED FOR THE WHOLE CHROMOSOME):
      LOH_hetSNP_number=floor(sum(LOHall$diff)/SNP_interval)
      print(paste("No. of Het SNPs to be added to LOH regions:",LOH_hetSNP_number))
    }
    # reconstruct allele counts for the LOH region based on actual depth for all to be perfect heterozygotes - allele counts remain as integers
    #m[cbind(1:nrow(m),2+m$a0)]=ifelse(m$depth%%2==0,m$depth/2,floor(m$depth/2))
    #m[cbind(1:nrow(m),2+m$a1)]=ifelse(m$depth%%2==0,m$depth/2,ceiling(m$depth/2))
    # the above two lines results in 4x number of Het SNPs compared with the rest of the chromosome -> need to have proportional number of Het SNPs in the LOH region
    ####
    lohs=data.frame() ####################################### get all non_LOH regions####
    for (j in 1:nrow(LOHall)){
      loh=ac[which(ac$position>=LOHall$start.pos[j] & ac$position<=LOHall$end.pos[j]),]
      m=merge(loh,al,"position")
      if (nrow(m)==nrow(loh)){
        print("merge OK")
      } else {print("ERROR - merge not OK")}
      # RE-reconstruct allele counts for LOH region
      hetSNP_number=LOHall$diff[j]/SNP_interval
      if (nrow(m)>hetSNP_number){
        print("more rows in LOH region than Het SNP number")
        for (k in 1:nrow(m)){
          if (k %% floor(nrow(m)/hetSNP_number)==0){
            m[cbind(k,2+m$a0[k])]=ifelse(m$depth[k]%%2==0,m$depth[k]/2,ceiling(m$depth[k]/2))
            m[cbind(k,2+m$a1[k])]=ifelse(m$depth[k]%%2==0,m$depth[k]/2,floor(m$depth[k]/2))
            print(k)
          }
        }
      } else {
        print("less rows in LOH region than Het SNP number - turning all into Heterozygotes")
        for (k in 1:nrow(m)){
          m[cbind(k,2+m$a0[k])]=ifelse(m$depth[k]%%2==0,m$depth[k]/2,ceiling(m$depth[k]/2))
          m[cbind(k,2+m$a1[k])]=ifelse(m$depth[k]%%2==0,m$depth[k]/2,floor(m$depth[k]/2))
          #print(k)
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
    } else {print("ERROR - missing SNPs - LOH and non-LOH regions not generated correctly; no AC file generated")}
  } else {
    ac_out=ac
    write.table(ac_out,paste0(NORMALNAME,"_alleleFrequencies_chr",i,".txt"),col.names=F,row.names=F,quote=F,sep="\t")
    print(paste("No changes made to the alleleCounter file - no LOH in chr",i))
  }
  #}
  ###
  ##
  #
  print(paste("STEP 2&3 - chr",i,"completed"))
}
