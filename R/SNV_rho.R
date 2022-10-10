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

#' Run ANNOVAR on somatic SNV loci to allow identification of germline variants
#'
#' The function assumes you have the latest version of ANNOVAR (Version 2019-10-24)
#' @param annovar_install_dir The full path to where table_annovar.pl can be found (if module load is used for annovar, use NULL).
#' @param humandb_dir The full path to where annovar human reference datasets can be found (i.e. gnomAD and other reference files)
#' @param input.vcf The name of the VCF file to be annotated (expected to be in working directory)
#' @param output.filename The name of the file to be used for the output (Default input.vcf)
#' @param genomebuild Genome build of the BAM file upon which the VCF file was obtained (options: "hg19" or "hg38") 
#' @param avsnp_version The version of dbSNP dataset of known SNPs (Default NULL)
#' @param exac_version The version of ExAC dataset for population frequencies of variants (Default NULL)
#' @param gnomad_version The version of gnomAD dataset for population frequencies of variants (Default gnomad211_genome)
#' @param clinvar_version The version of clinvar dataset for known pathogenic variants with clinical significance (Default NULL)
#' @author Naser Ansari-Pour (WIMM, Oxford)
#' @export

runANNOVAR = function(annovar_install_dir,humandb_dir,input.vcf,output.filename,genomebuild,avsnp_version,exac_version,gnomad_version,clinvar_version) {
  
  Operations=c(avsnp_version,exac_version,gnomad_version,clinvar_version)
  Noperation=length(Operations[!is.na(Operations)])
  
  cmd = paste(paste0(annovar_install_dir,"table_annovar.pl"),
              input.vcf,
              humandb_dir,
              "-buildver",
              genomebuild,
              "-out", output.filename,
              "-remove",
              "-protocol",
              paste("refGene",paste(Operations[which(!is.na(Operations))],collapse = ","),sep = ","),
              "-operation",
              paste("g",paste(rep("f",Noperation),collapse = ","),sep = ","),
              "-nastring .",
              "-vcfinput")

  print(cmd)
  
  system(cmd, wait=T)
  
  avinput_nrow = as.integer(system(paste0("wc -l ",input.vcf,".avinput | awk '{print $1}'"),intern = T))
  multianno_nrow = as.integer(system(paste0("wc -l ",input.vcf,".",genomebuild,"_multianno.txt | awk '{print $1}'"),intern = T))
  
  if ((avinput_nrow+1)==multianno_nrow){
    print("Successful annotation by ANNOVAR")
    VARtype <<-"vcf"
  } else {
    print("VCF file format error - restarting annotation from generated avinput")
    
    # clean up variants in avinput based on FILTER=="PASS" in the VCF file
    VCF=read.table(input.vcf,sep="\t",stringsAsFactors = F)
    avinput_file=read.table(paste0(input.vcf,".avinput"),sep="\t",stringsAsFactors = F)
    avinput_file=avinput_file[,1:5]
    if (nrow(avinput_file)==nrow(VCF)){
      avinput_pass=avinput_file[which(VCF$V7=="PASS"),]
      write.table(avinput_pass,paste0(input.vcf,".avinput"),col.names = F,row.names = F,quote = F,sep="\t")
    } else { stop("Clean-up of variants in avinput file failed - number of rows in avinput does not match variant count in VCF")}
    
    cmd_avinput = paste(paste0(annovar_install_dir,"table_annovar.pl"),
                paste0(input.vcf,".avinput"),
                humandb_dir,
                "-buildver",
                genomebuild,
                "-out", output.filename,
                "-remove",
                "-protocol",
                paste("refGene",paste(Operations[which(!is.na(Operations))],collapse = ","),sep = ","),
                "-operation",
                paste("g",paste(rep("f",Noperation),collapse = ","),sep = ","),
                "-nastring .")
    
    print(cmd_avinput)
    
    system(cmd_avinput, wait=T)
    
    multianno_nrow = as.integer(system(paste0("wc -l ",input.vcf,".",genomebuild,"_multianno.txt | awk '{print $1}'"),intern = T))
    
    if ((avinput_nrow+1)==multianno_nrow){
      print("Successful annotation by ANNOVAR from avinput file")
      VARtype <<- "avinput"
    } else {
      stop("Unknown VCF file format error - please regenerate VCF file")
    }
  }
}


#' Obtain an independent purity (rho) estimate based on the VAF distribution of somatic SNVs
#'
#' @param tumourname Identifier to be used for tumour output files (i.e. the tumour BAM file name without the '.bam' extension).
#' @param tumourbam Full path to the tumour BAM file
#' @param skip_balanced_region_finding Flag, set to TRUE if a prior balanced region bed file is available (Default FALSE)
#' @param chrom_coord Full path to the file with chromosome coordinates including start, end and left/right centromere positions
#' @param g1000alleles.prefix Prefix path to the 1000 Genomes SNP allele reference files
#' @param min_count Integer, minimum depth required for a SNP to be included in identifying hetSNPs (required, Default 10)
#' @param kmin The min number of SNPs to support a segment in PCF of 1000G hetSNP BAF values (Default 50)
#' @param gamma The PCF gamma value for segmentation of 1000G hetSNP BAF values (Default 150).
#' @param baf_tolerance The maximum deviation allowed for segmented BAF from 0.5 which is the expected mean for balanced regions (Default 0.02)
#' @param chrom_names A vector containing the names of chromosomes to be included (Default 1:22 or chr1,...,chr22)
#' @param skip_run_annovar Flag, set to TRUE if a prior somatic SNV file is already available (Default FALSE, a tumour-only VCF file is required in working directory)
#' @param genomebuild Genome build of the BAM file upon which the VCF file was obtained (options: "hg19" or "hg38") 
#' @param VCFprefix The character string that may have been assigned to the VCF file preceding the tumourname (Default=NULL)
#' @param VCFsuffix The character string that may have been assigned to the VCF file after the tumourname (Default=".vcf")
#' @param annovar_install_dir The full path to where table_annovar.pl can be found (if module load is used on the HPC cluster for annovar, use NULL).
#' @param humandb_dir The full path to where annovar human reference datasets can be found (i.e. gnomAD and other reference files)
#' @param maxAF The maximum gnomAD population allele frequency (AF) allowed for selecting somatic SNVs and removing SNPs (a value of 0 will remove all variants observed in gnomAD)  
#' @param avsnp_version The version of dbSNP dataset of known SNPs (Default NULL)
#' @param exac_version The version of ExAC dataset for population frequencies of variants (Default NULL)
#' @param gnomad_version The version of gnomAD dataset for population frequencies of variants (Default gnomad211_genome)
#' @param clinvar_version The version of clinvar dataset for known pathogenic variants with clinical significance (Default NULL)
#' @param VARprefix The character string that may have been assigned to the alleleCounter output file preceding the samplename (Default=NULL)
#' @param VARsuffix The character string that may have been assigned to the alleleCounter output file preceding the samplename (Default=NULL)
#' @param min_base_qual Minimum base quality required for a read to be counted when allele counting (Default: 20)
#' @param min_map_qual Minimum mapping quality required for a read to be counted when allele counting (Default: 35)
#' @param allelecounter_exe The full path to where the alleleCounter executable can be found (optional, default points to $PATH).
#' @param min_snv_depth Integer, minimum depth required for an SNV to be included in calculating rho (required, Default 10)
#' @param peak_threshold The minimum density required for a peak to be retained in the SNV VAF density distribution (Default 0.02)
#' @author Naser Ansari-Pour (WIMM, Oxford)
#' @export

getSNVrho=function(tumourname,tumourbam,skip_balanced_region_finding=FALSE,chrom_coord,g1000alleles.prefix,min_count=10,kmin=50,gamma=150,
                     baf_tolerance=0.02,chrom_names=1:22,skip_run_annovar=FALSE,genomebuild,VCFprefix,VCFsuffix,annovar_install_dir,
                     humandb_dir,maxAF,avsnp_version,exac_version,gnomad_version,clinvar_version,VARprefix=NULL,VARsuffix=NULL,
                     min_base_qual,min_map_qual,allelecounter_exe,min_snv_depth,peak_threshold=0.02){

  if (!skip_balanced_region_finding){
colClasses=c(chr="numeric",start="numeric",cen.left.base="numeric",cen.right.base="numeric",end="numeric")
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
  pcf_balanced=pcf_output[which(pcf_output$mean>(0.5-baf_tolerance) & pcf_output$mean<(0.5+baf_tolerance)),]
  balanced_regions=rbind(balanced_regions,pcf_balanced)
}
balanced_regions$length=balanced_regions$end.pos-balanced_regions$start.pos
print(paste("Proportion of genome in balanced state :",sum(balanced_regions$length)/sum(chr_loc$length)))
write.table(balanced_regions,paste0(tumourname,"_balanced_regions_for_snv_rho.txt"),col.names = T,row.names = F,quote = F,sep="\t")
}

  if (!skip_run_annovar){
runANNOVAR(annovar_install_dir=annovar_install_dir,
           humandb_dir=humandb_dir,
           input.vcf=paste0(VCFprefix,tumourname,VCFsuffix), 
           output.filename=paste0(VCFprefix,tumourname,VCFsuffix), 
           genomebuild=genomebuild,
           avsnp_version=avsnp_version,
           exac_version=exac_version,
           gnomad_version=gnomad_version,
           clinvar_version=clinvar_version)
    
# get variants from the variant file (VCF or ANNOVAR multianno.txt file)
  
  if (VARtype=="vcf"){
  VAR=data.frame(data.table::fread(paste0(VCFprefix,tumourname,VCFsuffix,".",genomebuild,"_multianno.txt"),stringsAsFactors=F)) # variant calls
  VARpass=VAR[VAR$Otherinfo10=="PASS",]
  } else if (VARtype=="avinput"){
    VARpass=data.frame(data.table::fread(paste0(VCFprefix,tumourname,VCFsuffix,".",genomebuild,"_multianno.txt"),stringsAsFactors=F)) # variant calls # already cleaned up in runANNOVAR
  }
  # remove germline variants with maxAF cutoff
  VARpass$AF[VARpass$AF=="."]=0
  VARpass$AF=as.numeric(VARpass$AF)
  if (maxAF>0){
    VARpass=VARpass[which(VARpass$AF<maxAF),]
  } else if (maxAF==0){
    VARpass=VARpass[which(VARpass$AF==maxAF),]
  }
  } else if (file.exists(paste0(VCFprefix,tumourname,VCFsuffix,".",genomebuild,"_multianno.txt"))) {
    VAR=data.frame(data.table::fread(paste0(VCFprefix,tumourname,VCFsuffix,".",genomebuild,"_multianno.txt"),stringsAsFactors=F)) # variant calls
    VARpass=VAR[VAR$Otherinfo10=="PASS",]
    VARpass$AF[VARpass$AF=="."]=0
    VARpass$AF=as.numeric(VARpass$AF)
    if (maxAF>0){
      VARpass=VARpass[which(VARpass$AF<maxAF),]
    } else if (maxAF==0){
      VARpass=VARpass[which(VARpass$AF==maxAF),]
    }
    print("Existing ANNOVAR multianno.txt file read")
  } else {
    VARpass=read.table(paste0(VARprefix,tumourname,VARsuffix),stringsAsFactors = F) # read in existing somatic variant file with five columns: chr, start, end, ref, alt
    print("Prior somatic variant file read")
  }
  # remove indels from variants
  VARpass$var=paste0(VARpass$Ref,VARpass$Alt)
  VARpass=VARpass[which(nchar(VARpass$var)==2 & !grepl("-",VARpass$var)),]
  
  # prepare input for allelecounting
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
