callChrXsubclones <- function(TUMOURNAME,X_GAMMA,X_KMIN,GENOMEBUILD,AR){
print(TUMOURNAME)
d=readLines(paste0(TUMOURNAME,"_mutantLogR_gcCorrected.tab"))
x=d[grepl("X",d)]
s=unlist(strsplit(x,split = "\t"))
X=matrix(s,nrow=length(x),byrow=T)
df=data.frame(chr="X",pos=as.numeric(X[-1,2]),logR=as.numeric(X[-1,3]))
df=df[which(df$pos>2600000 & df$pos<156000000),]
print(paste("dimension of chrX nonPAR logR file =",dim(df)))
PCF=pcf(df,gamma=X_GAMMA,kmin=X_KMIN)
write.table(PCF,paste0(TUMOURNAME,"_PCF_gamma_",X_GAMMA,"_chrX.txt"),col.names=T,row.names=F,quote=F,sep="\t")
print("PCF 1000 done")

if (GENOMEBUILD=="hg19"){
x_centromere=c(58632012,61632012) # hg19
ar=data.frame(startpos=66763874,endpos=66950461)
} else {
  x_centromere=c(58605580,62412542) #hg38
  ar=data.frame(startpos=67544021,endpos=67730619)
}

# INPUT for copy number inference
SAMPLEsegs=data.frame(PCF,stringsAsFactors=F)
pupl=read.table(paste0(TUMOURNAME,"_cellularity_ploidy.txt"),header=T,stringsAsFactors=F)
SAMPLEpurity=pupl$cellularity
#SAMPLEploidy=round(pupl$ploidy/2)*2
SAMPLEn=pupl$ploidy
print(paste(SAMPLEpurity,SAMPLEn))

# Estimating LogR deviation in diploid and gained regions (AUTOSOMAL)
BB=read.table(paste0(TUMOURNAME,"_subclones.txt"),header=T,stringsAsFactors = F)

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
print("likely WGD sample or no clonal LOH event")
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
      if (nrow(BBloh)>0){
        seg$clonal=ifelse(round(abs(explogrLoss-seg$mean),digits=2)<round(abs(sd(BBloh$LogR)/explogrLoss),digits=2),"yes","no")
      } else if (nrow(BBloh)==0){
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
    print(j)
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
      print(length(SPLIT[[j]]))
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
  outputDF$rank=NULL

print(paste("Number of rows merged =",nrow(SUBCLONESout)-nrow(outputDF)))
write.table(outputDF,paste0(TUMOURNAME,"_chrX_subclones.txt"),col.names = T,row.names = F,quote = F,sep="\t")

# PLOT
outputDF$diff=outputDF$endpos-outputDF$startpos
PGAclonal=sum(outputDF[which(outputDF$clonal=="yes"),]$diff)/sum(outputDF[which(!is.na(outputDF$clonal)),]$diff)


plot_BB=ggplot()+geom_hline(yintercept = 0:ceiling(max(outputDF$subclonalCN)),linetype="longdash",col="grey",size=0.2)+
  geom_rect(data=outputDF,aes(xmin=startpos,xmax=endpos,ymin=subclonalCN-0.02,ymax=subclonalCN+0.02))+
  geom_vline(xintercept = x_centromere,linetype="longdash",col="green")+
  #geom_hline(yintercept = nonpar,linetype="dotted",col="blue")+
  ylim(-0.2,ceiling(max(outputDF$subclonalCN))+0.2)+labs(x="ChrX coordinate (bp)",y="Average Ploidy")+
  theme(plot.title = element_text(hjust = 0.5,size=12),panel.background = element_blank())+
  ggtitle(paste0(TUMOURNAME," , PLOIDY: ",round(SAMPLEn,digits = 3)," , PURITY: ",round(SAMPLEpurity*100,digits = 0),
                 "%, PGAclonal: ",round(PGAclonal*100,digits = 1),"%"))

# ANDROGEN RECEPTOR LOCUS
if (AR){
  setDT(ar)
  setkey(ar,"startpos","endpos")
  setDT(outputDF)
  setkey(outputDF,"startpos","endpos")
  segAR=foverlaps(ar,outputDF,type="any",nomatch = 0)
  segAR$subclonalCN=(segAR$nMaj1+segAR$nMin1)*segAR$frac1+(segAR$nMaj2+segAR$nMin2)*segAR$frac2
  plot_BB=plot_BB+geom_rect(data=segAR,aes(xmin=startpos,xmax=endpos,ymin=subclonalCN-0.02,ymax=subclonalCN+0.02),fill="red")
}

pdf(paste0(TUMOURNAME,"_chrX_average_ploidy.pdf"))
print(plot_BB)
dev.off()

}
