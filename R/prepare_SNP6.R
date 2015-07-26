#' Adapted code from ASCAT to load in SNP6 data for plotting
#' noRD
ascat.loadData = function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, Germline_BAF_file = NULL, chrs = c(1:22,"X","Y"), gender = NULL, sexchromosomes = c("X","Y")) {
  
  # read in SNP array data files
  print.noquote("Reading Tumor LogR data...")
  Tumor_LogR <- read.table(Tumor_LogR_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
  print.noquote("Reading Tumor BAF data...")
  Tumor_BAF <- read.table(Tumor_BAF_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
  
  #infinite values are a problem - change those
  Tumor_LogR[Tumor_LogR==-Inf]=NA
  Tumor_LogR[Tumor_LogR==Inf]=NA
  
  Germline_LogR = NULL
  Germline_BAF = NULL
  if(!is.null(Germline_LogR_file)) {
    print.noquote("Reading Germline LogR data...")
    Germline_LogR <- read.table(Germline_LogR_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
    print.noquote("Reading Germline BAF data...")
    Germline_BAF <- read.table(Germline_BAF_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
    
    #infinite values are a problem - change those
    Germline_LogR[Germline_LogR==-Inf]=NA
    Germline_LogR[Germline_LogR==Inf]=NA
  }
  
  # make SNPpos vector that contains genomic position for all SNPs and remove all data not on chromosome 1-22,X,Y (or whatever is given in the input value of chrs)
  print.noquote("Registering SNP locations...")
  SNPpos <- Tumor_LogR[,1:2]
  SNPpos = SNPpos[SNPpos[,1]%in%chrs,]
  
  # if some chromosomes have no data, just remove them
  chrs = intersect(chrs,unique(SNPpos[,1]))
  
  Tumor_LogR = Tumor_LogR[,c(-1,-2),drop=F]
  Tumor_BAF = Tumor_BAF[,c(-1,-2),drop=F]
  # make sure it is all converted to numerical values
  for (cc in 1:dim(Tumor_LogR)[2]) {
    Tumor_LogR[,cc]=as.numeric(as.vector(Tumor_LogR[,cc]))
    Tumor_BAF[,cc]=as.numeric(as.vector(Tumor_BAF[,cc]))
  }
  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[,c(-1,-2),drop=F]
    Germline_BAF = Germline_BAF[,c(-1,-2),drop=F]
    for (cc in 1:dim(Germline_LogR)[2]) {
      Germline_LogR[,cc]=as.numeric(as.vector(Germline_LogR[,cc]))
      Germline_BAF[,cc]=as.numeric(as.vector(Germline_BAF[,cc]))
    }
  }
  
  # sort all data by genomic position
  last = 0;
  ch = list();
  SNPorder = vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke = SNPpos[SNPpos[,1]==chrs[i],]
    chrpos = chrke[,2]
    names(chrpos) = rownames(chrke)
    chrpos = sort(chrpos)
    ch[[i]] = (last+1):(last+length(chrpos))  
    SNPorder[ch[[i]]] = names(chrpos)
    last = last+length(chrpos)
  }
  SNPpos = SNPpos[SNPorder,]
  Tumor_LogR=Tumor_LogR[SNPorder,,drop=F]
  Tumor_BAF=Tumor_BAF[SNPorder,,drop=F]
  
  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[SNPorder,,drop=F]
    Germline_BAF = Germline_BAF[SNPorder,,drop=F]
  }
  
  # split the genome into distinct parts to be used for segmentation (e.g. chromosome arms, parts of genome between gaps in array design)
  print.noquote("Splitting genome in distinct chunks...")
  chr = split_genome(SNPpos)
  
  if (is.null(gender)) {
    gender = rep("XX",dim(Tumor_LogR)[2])
  }
  return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF, 
              Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL, 
              Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF, 
              SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs, 
              samples = colnames(Tumor_LogR), gender = gender, 
              sexchromosomes = sexchromosomes,
              failedarrays = NULL))
}


#' Parse the reference info file
#' @noRD
parseSNP6refFile = function(snp6_reference_info_file) {
  return(read.table(snp6_reference_info_file, header=T, stringsAsFactors=F))
}

#' Morph cel files into BAF and LogR
#'
cel2baf.logr = function(normal_cel_file, tumour_cel_file, output_file, snp6_reference_info_file, apt.probeset.genotype.exe="apt-probeset-genotype", apt.probeset.summarize.exe="apt-probeset-summarize", norm.geno.clust.exe="normalize_affy_geno_cluster.pl") {
  #PROBSET_GENOT = "~pvl/PennCNV/apt-1.12.0-20091012-i386-intel-linux/bin/apt-probeset-genotype" 
  #GW_SNP6 = "~pvl/PennCNV/gw6/lib/GenomeWideSNP_6.cdf"
  #SNP6_BIRDSEED_MODELS = "~pvl/PennCNV/gw6/lib/GenomeWideSNP_6.birdseed.models"
  #SNP6_SPECIALSNPS = "~pvl/PennCNV/gw6/lib/GenomeWideSNP_6.specialSNPs"
  
  #PROBSET_SUMMA = "~pvl/PennCNV/apt-1.12.0-20091012-i386-intel-linux/bin/apt-probeset-summarize"
  #QUANT_NORM_TARGET = "~pvl/PennCNV/gw6/lib/hapmap.quant-norm.normalization-target.txt"
  
  #NORM_GENO_CLUST = "~pvl/PennCNV/gw6/bin/normalize_affy_geno_cluster.pl"
  #UNM_NORMALS = "/nfs/team78pc2/pvl/normals/selected285/apt/gw6.genocluster"
  #LOCFILE = "~pvl/PennCNV/gw6/lib/affygw6.hg18.pfb"
  
  # Unpack pointers to reference files required during this step
  ref.files = parseSNP6refFile(snp6_reference_info_file)
  GW_SNP6 = ref.files[ref.files$variable == "GW_SNP6",]$reference_file
  SNP6_BIRDSEED_MODELS = ref.files[ref.files$variable == "SNP6_BIRDSEED_MODELS",]$reference_file
  SNP6_SPECIALSNPS = ref.files[ref.files$variable == "SNP6_SPECIALSNPS",]$reference_file
  QUANT_NORM_TARGET = ref.files[ref.files$variable == "QUANT_NORM_TARGET",]$reference_file
  LOCFILE = ref.files[ref.files$variable == "LOCFILE",]$reference_file
  UNM_NORMALS = ref.files[ref.files$variable == "UNM_NORMALS",]$reference_file
  
  # Unpack the normal cel file
  cmd = paste(apt.probeset.genotype.exe, "-c", GW_SNP6, "-a birdseed", "--read-models-birdseed", SNP6_BIRDSEED_MODELS, "--special-snps", SNP6_SPECIALSNPS, "--cels", normal_cel_file)
  print(cmd)
  system(cmd, wait=T)
  # Unpack the tumour cel file
  cmd = paste(apt.probeset.summarize.exe, "--cdf-file", GW_SNP6, "--analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true", "--target-sketch", QUANT_NORM_TARGET, normal_cel_file, tumour_cel_file)
  print(cmd)
  system(cmd, wait=T)
  # Construct the LogR and BAF and push that to 
  #cmd = paste(norm.geno.clust.exe, UNM_NORMALS, "quant-norm.pm-only.med-polish.expr.summary.txt", "-locfile", LOCFILE, "-out", paste(output_file_prefix,"_lrr_baf.txt", sep=""))
  cmd = paste(norm.geno.clust.exe, UNM_NORMALS, "quant-norm.pm-only.med-polish.expr.summary.txt", "-locfile", LOCFILE, "-out", output_file)
  print(cmd)
  system(cmd, wait=T)
}

#' Correct the BAF and LogR estimates for GC content
#'
gc.correct = function(samplename, infile.logr.baf, outfile.tumor.LogR, outfile.tumor.BAF, outfile.normal.LogR, outfile.normal.BAF, outfile.probeBAF, snp6_reference_info_file, chr_names, birdseed_report_file="birdseed.report.txt") {
  # Read in needed reference files
  ref.files = parseSNP6refFile(snp6_reference_info_file)
  SNP_POS_REF = ref.files[ref.files$variable == "SNP_POS",]$reference_file
  GC_SNP6 = ref.files[ref.files$variable == "GC_SNP6",]$reference_file
  
  #lrrbaf = read.table("lrr_baf.txt", header = T, sep = "\t", row.names=1)
  #lrrbaf = read.table(paste(lrr_baf_prefix, "_lrr_baf.txt", sep=""), header=T, sep="\t", row.names=1, stringsAsFactors=F)
  lrrbaf = read.table(infile.logr.baf, header=T, sep="\t", row.names=1, stringsAsFactors=F)
  
  #SNPpos = read.table("SNPpos.txt",header=T,sep="\t",row.names=1)
  #SNPpos = read.table(snp_pos_file, header=T, sep="\t", row.names=1, stringsAsFactors=F)
  SNPpos = read.table(SNP_POS_REF, header=T, sep="\t", row.names=1, stringsAsFactors=F)
  
  #sample = sub(".normal.CEL.Log.R.Ratio","",colnames(lrrbaf)[3])
  
  Tumor_LogR = lrrbaf[rownames(SNPpos), 5, drop=F]
  colnames(Tumor_LogR) = samplename
  
  Tumor_BAF = lrrbaf[rownames(SNPpos), 6, drop=F]
  colnames(Tumor_BAF) = samplename
  
  Normal_LogR = lrrbaf[rownames(SNPpos), 3, drop=F]
  colnames(Normal_LogR) = samplename
  
  Normal_BAF = lrrbaf[rownames(SNPpos), 4, drop=F]
  colnames(Normal_BAF) = samplename
  
  #replace 2's by NA
  Tumor_BAF[Tumor_BAF==2]=NA
  Normal_BAF[Normal_BAF==2]=NA
  
  # Tumor_LogR: correct difference between copy number only probes and other probes
  CNprobes = substring(rownames(SNPpos),1,2)=="CN"
  
  Tumor_LogR[CNprobes,1] = Tumor_LogR[CNprobes,1]-mean(Tumor_LogR[CNprobes,1],na.rm=T)
  Tumor_LogR[!CNprobes,1] = Tumor_LogR[!CNprobes,1]-mean(Tumor_LogR[!CNprobes,1],na.rm=T)  
  
  Normal_LogR[CNprobes,1] = Normal_LogR[CNprobes,1]-mean(Normal_LogR[CNprobes,1],na.rm=T)
  Normal_LogR[!CNprobes,1] = Normal_LogR[!CNprobes,1]-mean(Normal_LogR[!CNprobes,1],na.rm=T)  
  
  # limit the number of digits:
  Tumor_LogR = round(Tumor_LogR,4)
  Normal_LogR = round(Normal_LogR,4)
  
  #write.table(cbind(SNPpos,Tumor_BAF), paste(output_file_prefix, ".tumour.BAF.txt", sep=""), sep="\t", row.names=T, col.names=NA, quote=F)
  write.table(cbind(SNPpos,Tumor_BAF), paste(outfile.tumor.BAF, "_noGCcorr.txt", sep=""), sep="\t", row.names=T, col.names=NA, quote=F)
  #write.table(cbind(SNPpos,Normal_BAF), paste(output_file_prefix, ".normal.BAF.txt", sep=""), sep="\t", row.names=T, col.names=NA, quote=F)
  write.table(cbind(SNPpos,Normal_BAF), paste(outfile.normal.BAF, "_noGCcorr.txt", sep=""), sep="\t", row.names=T, col.names=NA, quote=F)
  
  # read into ASCAT and make GC corrected input:
  #write.table(cbind(SNPpos,Tumor_LogR), paste(output_file_prefix, ".tumour.LogR.txt", sep=""), sep="\t", row.names=T, col.names=NA, quote=F)
  write.table(cbind(SNPpos,Tumor_LogR), paste(outfile.tumor.LogR, "_noGCcorr.txt", sep=""), sep="\t", row.names=T, col.names=NA, quote=F)
  #write.table(cbind(SNPpos,Normal_LogR), paste(output_file_prefix, ".normal.LogR.txt", sep=""), sep="\t", row.names=T, col.names=NA, quote=F)
  write.table(cbind(SNPpos,Normal_LogR), paste(outfile.normal.LogR, "_noGCcorr.txt", sep=""), sep="\t", row.names=T, col.names=NA, quote=F)
  
  # ======================================= above previous prepareGCcorrect, below runGCcorrect ==============================================
  
  # TODO: This must be a dapted to not hardcode the chromosome names
  gender <- read.table(birdseed_report_file, sep="\t", skip=66, header=T)
  sex <- as.vector(gender[,"computed_gender"])
  sex[sex == "female"] <- "XX"
  sex[sex == "male"] <- "XY"
  sex[sex == "unknown"] <- NA
  
  ascat.bc <- ascat.loadData(paste(outfile.tumor.LogR, "_noGCcorr.txt", sep=""), paste(outfile.tumor.BAF, "_noGCcorr.txt", sep=""),paste(outfile.normal.LogR, "_noGCcorr.txt", sep=""), paste(outfile.normal.BAF, "_noGCcorr.txt", sep=""), chrs=chr_names, gender=sex)
  ascat.plotRawData(ascat.bc)
  ascat.bc <- ascat.GCcorrect(ascat.bc, GC_SNP6)

  # Make sure the right column names are added here, because these are expected by fitcopynumber
  colnames(ascat.bc$SNPpos) = c("Chromosome", "Position")


# Determine SNPs with BAF between 0.3-0.7 from normal => these are supposed to be heterozygous
  is.het = (ascat.bc$Germline_BAF >= 0.3 & ascat.bc$Germline_BAF <= 0.7)
  dat = cbind(ascat.bc$SNPpos, round(ascat.bc$Germline_LogR, 4))
  dat = dat[which(is.het),]
  colnames(dat) = c("Chromosome", "Position", samplename)
  write.table(dat, file=outfile.normal.LogR, row.names=T, quote=F, sep="\t")

  dat = cbind(ascat.bc$SNPpos, round(ascat.bc$Germline_BAF, 4))
  dat = dat[which(is.het),]
  colnames(dat) = c("Chromosome", "Position", samplename)
  write.table(dat, file=outfile.normal.BAF, row.names=T, quote=F, sep="\t")

  # Save the probe ids plus their BAF for only the germline heterozygous mutations
  select = !is.na(ascat.bc$Tumor_BAF)
  dat = cbind(row.names(ascat.bc$SNPpos), ascat.bc$Tumor_BAF)
  dat = dat[which(select & is.het),]
  write.table(dat, file=outfile.probeBAF, row.names=F, quote=F, col.names=F, sep="\t")

  dat = cbind(ascat.bc$SNPpos, round(ascat.bc$Tumor_BAF, 4))
  dat = dat[which(select),]
  colnames(dat) = c("Chromosome", "Position", samplename)
  write.table(dat, file=outfile.tumor.BAF, row.names=T, quote=F, sep="\t")

  # Save tumour BAF and LogR directly. Include homozygous SNPs here.
  select = !is.na(ascat.bc$Tumor_LogR)
  dat = cbind(ascat.bc$SNPpos, round(ascat.bc$Tumor_LogR, 4))
  dat = dat[which(select),]
  colnames(dat) = c("Chromosome", "Position", samplename)
  write.table(dat, file=outfile.tumor.LogR, row.names=T, quote=F, sep="\t")
}

#' Prepares data for impute
#'
generate.impute.input.snp6 = function(infile.probeBAF, outFileStart, chrom, chr_names, problemLociFile, snp6_reference_info_file, imputeinfofile, is.male, heterozygousFilter="none") {
  # Obtain pointer to SNP6 specific reference file
  ref.files = parseSNP6refFile(snp6_reference_info_file)
  ANNO_FILE = ref.files[ref.files$variable == "ANNO_FILE",]$reference_file
  
  # Read in the 1000 genomes reference file paths for the specified chrom
  impute.info = parse.imputeinfofile(imputeinfofile, is.male, chrom=chrom)
  #chr_names = unique(impute.info$chrom)
  
  #print(paste("GenerateImputeInput is.male? ", is.male,sep=""))
  #print(paste("GenerateImputeInput #impute files? ", nrow(impute.info),sep=""))
  
  # Read in the known SNP locations from the 1000 genomes reference files
  known_SNPs = read.table(impute.info$impute_legend[1], sep=" ", header=T)
  if(nrow(impute.info)>1){
    for(r in 2:nrow(impute.info)){
      known_SNPs = rbind(known_SNPs, read.table(impute.info$impute_legend[r], sep=" ", header=T))
    }
  }
  
  outfile=paste(outFileStart,chrom,".txt",sep="")
#   chr_names=c(1:22,"X")
#   if(chr==23){
#     #known_SNPs<-read.table(paste("/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL_1000G_phase1integrated_feb2012_impute/ALL_1000G_phase1integrated_feb2012_chrX_PAR1_impute.legend",sep=""),sep=" ",header=T)
#     known_SNPs<-read.table(known_snps_x_par1_file,sep=" ",header=T)
#     print(paste("#known SNPs PAR1=",nrow(known_SNPs),sep=""))
#     #known_SNPs<-rbind(known_SNPs,read.table(paste("/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL_1000G_phase1integrated_feb2012_impute/ALL_1000G_phase1integrated_feb2012_chrX_nonPAR_impute.legend",sep=""),sep=" ",header=T))
#     known_SNPs<-rbind(known_SNPs, read.table(known_snps_x_nonpar_file,sep=" ",header=T))
#     print(paste("#known SNPs PAR1 + nonPAR=",nrow(known_SNPs),sep=""))
#     #known_SNPs<-rbind(known_SNPs,read.table(paste("/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL_1000G_phase1integrated_feb2012_impute/ALL_1000G_phase1integrated_feb2012_chrX_PAR2_impute.legend",sep=""),sep=" ",header=T))
#     known_SNPs<-rbind(known_SNPs, read.table(known_snps_x_par2_file,sep=" ",header=T))
#     print(paste("#known SNPs PAR1 + nonPAR + PAR2=",nrow(known_SNPs),sep=""))
#   }else{
#     #known_SNPs<-read.table(paste("/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL_1000G_phase1integrated_feb2012_impute/ALL_1000G_phase1integrated_feb2012_chr",chr_names[chr],"_impute.legend",sep=""),sep=" ",header=T)
#     known_SNPs<-read.table(known_snps_autosomes_file,sep=" ",header=T)
#   }
  
  #280412
  #known_SNPs = known_SNPs[known_SNPs[,5]=="SNP",]
  
  #if(chr==23){
  #	known_SNPs<-read.table(paste("/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL.chrX.phase1.legend",sep=""),sep=" ",header=T,stringsAsFactors=F)
  #}else{
  #	known_SNPs<-read.table(paste("/lustre/scratch110/sanger/dw9/haplotype_pipeline/impute/ALL_1000G_phase1interim_jun2011_impute/ALL_1000G_phase1interim_jun2011_chr",chr,"_impute.legend",sep=""),sep=" ",header=T,stringsAsFactors=F)
  #}
  
  #make sure all bases are repesented as factors, in the correct order
  #known_SNPs$allele0=factor(known_SNPs$allele0,levels=c("A","C","G","T"))
  #known_SNPs$allele1=factor(known_SNPs$allele1,levels=c("A","C","G","T"))
  known_SNPs[,3]=factor(known_SNPs[,3],levels=c("A","C","G","T"))
  known_SNPs[,4]=factor(known_SNPs[,4],levels=c("A","C","G","T"))
#   chr_name=chrom
#   if(chr_name=="23")
#   {
#     chr_name="X"
#   }
  print(head(known_SNPs))
  print(dim(known_SNPs))
  chr_name = chr_names[chrom]
  
  #filter out bad SNPs (streaks in BAF)
  if((problemLociFile !="NA") & (!is.na(problemLociFile)))
  {
    problemSNPs=read.table(problemLociFile,header=T,sep="\t")
    problemSNPs=problemSNPs$Pos[problemSNPs$Chr==chr_name]
    badIndices=match(known_SNPs[,2],problemSNPs)
    known_SNPs = known_SNPs[is.na(badIndices),]
    print(paste("badIndices lengths=",length(badIndices),",",sum(is.na(badIndices)),sep=""))
  }
  
  #knownSNP6data=read.csv("/lustre/scratch110/sanger/dw9/haplotype_pipeline/GenomeWideSNP_6.na32.annot.csv",comment.char="#",header=T,row.names=NULL,stringsAsFactors=F)
  knownSNP6data=read.csv(ANNO_FILE,comment.char="#",header=T,row.names=NULL,stringsAsFactors=F)
  knownSNP6data=knownSNP6data[knownSNP6data$Chromosome==chr_name,]
  print(paste("first column=",names(knownSNP6data)[1],sep=""))
  print(paste("first known datum=",knownSNP6data[1,1],sep=""))
  
  #adjust for strand
  knownSNP6data$Allele.A[knownSNP6data$Strand=="-" & knownSNP6data$Allele.A=="A"]="X"
  knownSNP6data$Allele.A[knownSNP6data$Strand=="-" & knownSNP6data$Allele.A=="C"]="Y"
  knownSNP6data$Allele.A[knownSNP6data$Strand=="-" & knownSNP6data$Allele.A=="G"]="Z"
  knownSNP6data$Allele.A[knownSNP6data$Strand=="-" & knownSNP6data$Allele.A=="T"]="A"
  knownSNP6data$Allele.A[knownSNP6data$Strand=="-" & knownSNP6data$Allele.A=="X"]="T"
  knownSNP6data$Allele.A[knownSNP6data$Strand=="-" & knownSNP6data$Allele.A=="Y"]="G"
  knownSNP6data$Allele.A[knownSNP6data$Strand=="-" & knownSNP6data$Allele.A=="Z"]="C"
  knownSNP6data$Allele.B[knownSNP6data$Strand=="-" & knownSNP6data$Allele.B=="A"]="X"
  knownSNP6data$Allele.B[knownSNP6data$Strand=="-" & knownSNP6data$Allele.B=="C"]="Y"
  knownSNP6data$Allele.B[knownSNP6data$Strand=="-" & knownSNP6data$Allele.B=="G"]="Z"
  knownSNP6data$Allele.B[knownSNP6data$Strand=="-" & knownSNP6data$Allele.B=="T"]="A"
  knownSNP6data$Allele.B[knownSNP6data$Strand=="-" & knownSNP6data$Allele.B=="X"]="T"
  knownSNP6data$Allele.B[knownSNP6data$Strand=="-" & knownSNP6data$Allele.B=="Y"]="G"
  knownSNP6data$Allele.B[knownSNP6data$Strand=="-" & knownSNP6data$Allele.B=="Z"]="C"
  
  #remove duplicates (variants on both strands)
  knownSNP6data = knownSNP6data[!duplicated(knownSNP6data$Physical.Position),]
  
  #make sure all bases are repesented as factors, in the correct order
  knownSNP6data$Allele.A = factor(knownSNP6data$Allele.A,levels=c("A","C","G","T"))
  knownSNP6data$Allele.B = factor(knownSNP6data$Allele.B,levels=c("A","C","G","T"))
  
  snp_data = read.table(infile.probeBAF,sep="\t",header=F,row.names=NULL,stringsAsFactors=F)
  print(paste("first datum=",snp_data[1,1],sep=""))
  
  indices = match(snp_data[,1],knownSNP6data$Probe.Set.ID)
  if(sum(!is.na(indices))==0){
    indices = match(snp_data[,1],knownSNP6data$dbSNP.RS.ID)
  }
  print(paste("found SNPs=",sum(!is.na(indices)),sep=""))
  print(paste("class=",class(knownSNP6data$Physical.Position),sep=""))
  matched.info = cbind(knownSNP6data[indices[!is.na(indices)],c("Physical.Position","Allele.A","Allele.B")],snp_data[!is.na(indices),2])
  #matched.info=cbind(as.numeric(knownSNP6data[indices[!is.na(indices)],"Physical.Position"]),knownSNP6data[indices[!is.na(indices)],c("Allele.A","Allele.B")],snp_data[!is.na(indices),2])
  print(paste("class2=",class(matched.info[,1]),sep=""))
  #make sure all bases are repesented as factors, in the correct order
  #matched.info$Allele.A=factor(matched.info$Allele.A,levels=c("A","C","G","T"))
  #matched.info$Allele.B=factor(matched.info$Allele.B,levels=c("A","C","G","T"))
  
  print(paste("first row of matched.info=",paste(matched.info[1,],sep=","),sep=""))
  print(paste("first Allele.A=",matched.info$Allele.A[1],sep=""))
  print(paste("first Allele.B=",matched.info$Allele.B[1],sep=""))
  
  print(paste("class 1a =",class(known_SNPs[,2]),sep=""))
  print(paste("class 2a =",class(as.numeric(known_SNPs[,2])),sep=""))
  
  indices2 = match(matched.info[,1],known_SNPs[,2])
  #indices2<-match(matched.info[,1],as.numeric(known_SNPs[,2]))
  #combined.info=cbind(matched.info[!is.na(indices2),],known_SNPs[indices2[!is.na(indices2)],])
  combined.info = cbind(matched.info[!is.na(indices2),],known_SNPs[indices2[!is.na(indices2)],1:4])
  print(paste("first row of combined.info=",paste(combined.info[1,],sep=","),sep=""))
  lev2 = levels(combined.info[,2])
  print(paste("levels[2]=",paste(lev2,sep=","),sep=""))
  lev3 = levels(combined.info[,3])
  print(paste("levels[3]=",paste(lev3,sep=","),sep=""))
  lev7 = levels(combined.info[,7])
  print(paste("levels[7]=",paste(lev7,sep=","),sep=""))
  lev8 = levels(combined.info[,8])
  print(paste("levels[8]=",paste(lev8,sep=","),sep=""))
  
  combined.info1 = combined.info[(combined.info[,2]==combined.info[,7] & combined.info[,3]==combined.info[,8]),]
  #alleles are reversed
  combined.info2 = cbind(combined.info[(combined.info[,2]==combined.info[,8] & combined.info[,3]==combined.info[,7]),1:3],1.0-combined.info[(combined.info[,2]==combined.info[,8] & combined.info[,3]==combined.info[,7]),4],combined.info[(combined.info[,2]==combined.info[,8] & combined.info[,3]==combined.info[,7]),5:8])
  names(combined.info2) = names(combined.info1)
  
  all.info = rbind(combined.info1,combined.info2)
  
  print(paste("norows all.info=",nrow(all.info),sep=""))
  
  #all.info=all.info[order(all.info[,1]),]
  all.info = all.info[order(as.numeric(all.info[,1])),]
  names(all.info)[4]="allele.frequency"
  write.csv(all.info, file=paste(outFileStart,chrom,"_withAlleleFreq.csv",sep=""), quote=F, row.names=F)
  
  #no.rows=1:sum(!is.na(indices))
  no.rows = nrow(all.info)
  snp.names = paste("snp",1:no.rows,sep="")
  #out.data<-cbind(snp.names,all.info[!is.na(indices),5:8],matrix(data=c(0,1,0),nrow=no.rows,ncol=3,byrow=T))
  out.data = cbind(snp.names,all.info[5:8], matrix(data=c(0,1,0), nrow=no.rows, ncol=3, byrow=T))
  write.table(out.data,file=outfile,row.names=F,col.names=F,quote=F)
  if(chrom==23){
    sample.g.file = paste(outFileStart,"sample_g.txt",sep="")
    sample_g_data = data.frame(ID_1=c(0,"INDIVI1"),ID_2=c(0,"INDIVI1"),missing=c(0,0),sex=c("D",2))
    write.table(sample_g_data, file=sample.g.file, row.names=F, col.names=T, quote=F)
  }
}

#' Infer the gender using the birdseed report file
infer_gender_birdseed = function(birdseed_report_file) {
  z = read.table(birdseed_report_file, header=T)
  return(as.character(z$em.cluster.chrX.het.contrast_gender))
}
