#' Match allele Counts (MaC) (source: https://github.com/nansari-pour/MaC)
#' 
#' Function to match the allelic labels with the nucleotide counts obtained by alleleCounter at a given single nucleotide locus 
#' This function will calculate the depth and variant/B-allele allele frequency (vaf/baf) for SNV/SNP
#' It requires the accompanying alleles file for the loci file used in alleleCounter which has four columns for
#' chromosome, position, reference nucleotide and alternative nucleotide of SNV/SNP and is assumed to have a header for these columns
#' @param samplename The sample name used in allelecounting - samplename is usually the the BAM ID (e.g. you'll have samplename.bam)
#' @param allelecount.prefix The character string that may have been assigned to the alleleCounter output file preceding the samplename (Default=NULL)
#' @param allelecount.suffix The character string that may have been assigned to the alleleCounter output file following the samplename (Default="_ac.txt")
#' @param allelesfile.prefix The character string that may have been assigned to the alleles file (cotaining ref and alt nucleotides of SNV/SNP) output file preceding the samplename (Default=NULL)
#' @param allelesfile.suffix The character string that may have been assigned to the alleles file (cotaining ref and alt nucleotides of SNV/SNP) output file following the samplename (Default="_alleles.txt")
#' @param output.suffix The character string to be added following the samplename to the output file name (Default="_mac.txt")
#' @param remove.na Should the variants with an NA vaf be removed i.e. SNV/SNP with depth=0? (Default=TRUE, use FALSE if output is required at same length as alleleCounter output)
#' @author Naser Ansari-Pour (WIMM, Oxford)
#' @export

MaC=function(samplename,allelecount.prefix,allelecount.suffix,allelesfile.prefix,allelesfile.suffix,output.suffix,remove.na=TRUE){
  ac=read.table(paste0(allelecount.prefix,samplename,allelecount.suffix),stringsAsFactors=F)
  notation=ac$V1[1]
  print(paste("Chr notation example :",notation))
  if (nchar(notation)>3){
    ac$V1=gsub("chr","",ac$V1)
    ac$V1[ac$V1=="X"]=23
    ac$V1=as.numeric(ac$V1)
    ac=ac[order(ac$V1,ac$V2),]
  } else {
    ac$V1[ac$V1=="X"]=23
    ac$V1=as.numeric(ac$V1)
    ac=ac[order(ac$V1,ac$V2),]
  }
  print(paste("No. of variants with counts =",nrow(ac)))
  al=read.table(paste0(allelesfile.prefix,samplename,allelesfile.suffix),header=T,stringsAsFactors=F)
  colnames(al)=c("chr","pos","ref","alt")
  if (nchar(notation)>3){
    al$chr=gsub("chr","",al$chr)
    al$chr[al$chr=="X"]=23
    al$chr=as.numeric(al$chr)
    al=al[order(al$chr,al$pos),]
  } else {
    al$chr[al$chr=="X"]=23
    al$chr=as.numeric(al$chr)
    al=al[order(al$chr,al$pos),]
  }
  if (nrow(al[which(al$ref=="-" | al$alt=="-"),])>0){
    stop("Indels (insertion/deletion variants) detected in the alleles file - MaC only works for single nucleotide substitutions (i.e. SNV/SNP)")
  }
  al$ref[al$ref=="A"]=1
  al$ref[al$ref=="C"]=2
  al$ref[al$ref=="G"]=3
  al$ref[al$ref=="T"]=4
  al$alt[al$alt=="A"]=1
  al$alt[al$alt=="C"]=2
  al$alt[al$alt=="G"]=3
  al$alt[al$alt=="T"]=4
  ref=as.numeric(al$ref)
  ref_df=data.frame(pos=1:nrow(al),ref=ref+2)
  REF=ac[cbind(ref_df$pos,ref_df$ref)]
  alt=as.numeric(al$alt)
  alt_df=data.frame(pos=1:nrow(al),alt=alt+2)
  ALT=ac[cbind(alt_df$pos,alt_df$alt)]
  mac=data.frame(chr=al$chr,pos=al$pos,ref=as.numeric(REF),alt=as.numeric(ALT))
  mac$depth=mac$ref+mac$alt
  mac$vaf=mac$alt/mac$depth
  if (remove.na){
  mac=mac[which(!is.na(mac$vaf)),]
  print(paste("No. of variants retained after removing vaf==NA is",nrow(mac)))
  }
  # Revert chromosome X label (23 -> X)
  mac$chr[mac$chr==23]="X"
  # If Chr notation has 'chr', add it back
  if (nchar(notation)>3){
    mac$chr=paste0("chr",mac$chr)
  }
  write.table(mac,paste0(samplename,output.suffix),col.names=T,row.names=F,quote=F,sep="\t")
}


#' Obtain allele counts for somatic SNV loci to allow purity estimation through external program alleleCounter
#'
#' @param bam.file The BAM file of the sample
#' @param output.file The filename where count output will be written to.
#' @param g1000.loci A file with 1000 Genomes SNP loci.
#' @param min.base.qual The minimum base quality required for a base be counted (optional, default=20).
#' @param min.map.qual The minimum mapping quality required for a read to be counted (optional, default=35).
#' @param allelecounter.exe The full path to where the alleleCounter executable can be found (optional, default points to $PATH).
#' @author sd11, Naser Ansari-Pour (WIMM, Oxford)
#' @export

getSNValleleCounts = function(bam.file, output.file, loci, min.base.qual=20, min.map.qual=35, allelecounter.exe="alleleCounter") {
  cmd = paste(allelecounter.exe,
              "-b", bam.file,
              "-l", loci,
              "-o", output.file,
              "-m", min.base.qual,
              "-q", min.map.qual)
  
  
  # alleleCount >= v4.0.0 is sped up considerably on 1000G loci when run in dense-snp mode
  counter_version = system(paste(allelecounter.exe, "--version"), intern = T)
  if (as.integer(substr(x = counter_version, start = 1, stop = 1)) >= 4)
    cmd = paste(cmd, "--dense-snps")
  print(cmd)
  
  system(cmd, wait=T)
}

#' Obtain an independent purity (rho) estimate based on the VAF distribution of somatic SNVs
#'
#' @param tumourbam The BAM file of the sample
#' @param To be filled (NAP)
#' @param 
#' @param 
#' @param 
#' @param allelecounter.exe The full path to where the alleleCounter executable can be found (optional, default points to $PATH).
#' @author Naser Ansari-Pour (WIMM, Oxford)

getSNVrho=function(tumourname,tumourbam,chrom_coord,min_count=10,kmin=50,gamma=150,
                     baf_tolerance=0.02,chrom_names,VARtype,VARprefix,VARsuffix,maxAF,
                     gnomAD=T,min_base_qual,min_map_qual,allelecounter.exe,min_snv_depth=10,peak_threshold=0.02){

chr_loc=read.table(chrom_coord,colClasses = colClasses,header=T,stringsAsFactors = F) # chrom_coord = full path to chromosome coordinates 
chr_loc$length=(chr_loc$cen.left.base-chr_loc$start)+(chr_loc$end-chr_loc$cen.right.base)

balanced_regions=data.frame()
for (chr in chrom_names){
  print(chr)
  # read in alleleCounter output for each chromosome
  ac=read.table(paste0(tumourname,"_alleleFrequencies_chr",chr,".txt"),stringsAsFactors = F)
  ac=ac[order(ac$V2),]
  # match allele counts with respective SNP alleles
  al=read.table(paste0(g1000alleles.prefix,chr,".txt"),header=T,stringsAsFactors = F)
  ref=al$a0
  ref_df=data.frame(pos=1:nrow(al),ref=ref+2)
  REF=ac[cbind(ref_df$pos,ref_df$ref)]
  alt=al$a1
  alt_df=data.frame(pos=1:nrow(al),alt=alt+2)
  ALT=ac[cbind(alt_df$pos,alt_df$alt)]
  mac=data.frame(ref=REF,alt=ALT)
  mac$depth=as.numeric(mac$ref)+as.numeric(mac$alt)
  mac$baf=as.numeric(mac$alt)/as.numeric(mac$depth)
  mac=cbind(al,mac)
  mac=mac[which(mac$depth>=min_count),]
  hetmac=mac[which(mac$baf>0.1 & mac$baf<0.9),]
  pcf_input=data.frame(chr=chr,hetmac[,c("position","baf")])
  pcf_output=copynumber::pcf(pcf_input,kmin = kmin,gamma=gamma)
  pcf_balanced=pcf_output[which(pcf_output$mean>0.5-baf_tolerance & pcf_output$mean<0.5+baf_tolerance),]
  balanced_regions=rbind(balanced_regions,pcf_balanced)
  balanced_regions$length=balanced_regions$end.pos-balanced_regions$start.pos
}
print(paste("Proportion of genome in balanced state :",sum(balanced_regions$length)/sum(chr_loc$length)))
write.table(balanced_regions,paste0(tumourname,"_balanced_regions_for_snv_rho.txt"),col.names = T,row.names = F,quote = F,sep="\t")

# get variants from the variant file (VCF or ANNOVAR multianno.txt file)
if (VARtype=="ANNOVAR"){
  VAR=data.frame(data.table::fread(paste0(VARprefix,tumourname,VARsuffix),stringsAsFactors=F)) # variant calls
  if (length(which(grepl("Otherinfo",names(VAR))))>0){
  VARpass=VAR[VAR$Otherinfo10=="PASS",]
  }
  if (gnomAD){
  VARpass$AF[VARpass$AF=="."]=0
  VARpass$AF=as.numeric(VARpass$AF)
  VARpass=VARpass[which(VARpass$AF<maxAF),]
  }
  # remove indels from variants
  VARpass$var=paste0(VARpass$Ref,VARpass$Alt)
  VARpass=VARpass[which(nchar(VARpass$var)==2 & !grepl("-",VARpass$var)),]
} else if (VARtype=="VCF"){
  VAR=data.frame(data.table::fread(paste0(VARprefix,tumourname,VARsuffix),sep="\t", header=TRUE, skip="#CHROM")) # variant calls
  VARpass=VAR[VAR$FILTER=="PASS",]
  VARpass$var=paste0(VARpass$REF,VARpass$ALT)
  VARpass=VARpass[which(nchar(VARpass$var)==2 & !grepl("-",VARpass$var)),]
} else {stop("Unknown variant file type")}

VARpass=VARpass[,c(1,2,4,5)]
names(VARpass)=c("chr","pos","ref","alt")
notation=VARpass$chr[1]
print(paste("Chr notation example :",notation))
if (nchar(notation)>3){
VARpass$chr=gsub("chr","",VARpass$chr)
}
VARpass=VARpass[which(!is.na(match(VARpass$chr,chrom_names))),]
VARpass$chr=as.numeric(VARpass$chr)
VARpass=VARpass[order(VARpass$chr,VARpass$pos),]

# overlap with balanced regions
balanced_regions=read.table(paste0(tumourname,"_balanced_regions_for_snv_rho.txt"),header=T,stringsAsFactors = F)
balanced_regions=data.frame(chr=balanced_regions$chrom,startpos=balanced_regions$start.pos,endpos=balanced_regions$end.pos,stringsAsFactors = F)
# select variants in balanced regions
VARpass=VARpass[,c("chr","pos","pos","ref","alt")]
names(VARpass)=c("chr","startpos","endpos","ref","alt")
data.table::setDT(VARpass)
data.table::setDT(balanced_regions)
data.table::setkey(balanced_regions,"chr","startpos","endpos")
data.table::setkey(VARpass,"chr","startpos","endpos")
VARpass_balanced=data.table::foverlaps(VARpass,balanced_regions,type="within",nomatch = 0)
VARfinal=VARpass_balanced[,c(1,5:7)]
names(VARfinal)=c("chr","pos","ref","alt")

if (nchar(notation)>3){
  VARfinal$chr=paste0("chr",VARfinal$chr)
}
print(head(VARfinal))
write.table(VARfinal,paste0(tumourname,"_alleles.txt"),col.names=T,row.names=F,quote=F,sep="\t")
write.table(VARfinal[,1:2],paste0(tumourname,"_loci.txt"),col.names=F,row.names=F,quote=F,sep="\t")

# Run alleleCounter
getSNValleleCounts(bam.file=tumourbam,
                output.file=paste0(tumourname,"_ac.txt"),
                loci=paste0(tumourname,"_loci.txt"),
                min.base.qual=min_base_qual,
                min.map.qual=min_map_qual,
                allelecounter.exe=allelecounter_exe)

# Run Match alleleCounts
MaC(samplename=tumourname,
    allelecount.prefix = NULL,
    allelecount.suffix = "_ac.txt",
    allelesfile.prefix = NULL,
    allelesfile.suffix = "_alleles.txt",
    output.suffix = "_mac.txt",
    remove.na = TRUE)

# Get density of VAF based on somatic SNVs in balanced regions
mac=read.table(paste0(tumourname,"_mac.txt"),header = T,stringsAsFactors = F)
VAF=mac[mac$depth>=min_snv_depth,]
x=density(VAF$vaf)
y=x$y
 
# peaks have to be at least peak_threshold% of the largest peak in the density distribution 
peak_density=NULL
for (i in 2:(length(y)-1)){
  if (y[i]>y[(i-1)] & y[i]>y[(i+1)]){
    peak_density=append(peak_density,y[i])
    print(paste0(i,"th = ",y[i]))
  }
}
peak_density=peak_density[which(peak_density>max(peak_density)*peak_threshold)]
peak_vaf=x$x[which(!is.na(match(x$y,peak_density)))]

D=ggplot(VAF,aes(vaf))+geom_density(col="green4",fill="green4")+geom_vline(xintercept = peak_vaf,col="red")+xlab("VAF")+xlim(0,1)
ggsave(paste0(tumourname,"_VAF_peaks.pdf"),device="pdf")
write.table(data.frame(vaf_peaks=round(peak_vaf,digits = 2),stringsAsFactors = F),paste0(tumourname,"_VAF_peaks.txt"),col.names = T,row.names = F,quote = F)
write.table(data.frame(snv_rho=round(2*max(peak_vaf),digits = 2),stringsAsFactors = F),paste0(tumourname,"_snv_rho.txt"),col.names = T,row.names = F,quote = F)
}

