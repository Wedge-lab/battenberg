#' Adapted code from ASCAT to load in SNP6 data for plotting
#' noRD
# ascat.loadData = function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, Germline_BAF_file = NULL, chrs = c(1:22,"X","Y"), gender = NULL, sexchromosomes = c("X","Y")) {
#   
#   # read in SNP array data files
#   print.noquote("Reading Tumor LogR data...")
#   Tumor_LogR <- read.table(Tumor_LogR_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
#   print.noquote("Reading Tumor BAF data...")
#   Tumor_BAF <- read.table(Tumor_BAF_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
#   
#   #infinite values are a problem - change those
#   Tumor_LogR[Tumor_LogR==-Inf]=NA
#   Tumor_LogR[Tumor_LogR==Inf]=NA
#   
#   Germline_LogR = NULL
#   Germline_BAF = NULL
#   if(!is.null(Germline_LogR_file)) {
#     print.noquote("Reading Germline LogR data...")
#     Germline_LogR <- read.table(Germline_LogR_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
#     print.noquote("Reading Germline BAF data...")
#     Germline_BAF <- read.table(Germline_BAF_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
#     
#     #infinite values are a problem - change those
#     Germline_LogR[Germline_LogR==-Inf]=NA
#     Germline_LogR[Germline_LogR==Inf]=NA
#   }
#   
#   # make SNPpos vector that contains genomic position for all SNPs and remove all data not on chromosome 1-22,X,Y (or whatever is given in the input value of chrs)
#   print.noquote("Registering SNP locations...")
#   SNPpos <- Tumor_LogR[,1:2]
#   SNPpos = SNPpos[SNPpos[,1]%in%chrs,]
#   
#   # if some chromosomes have no data, just remove them
#   chrs = intersect(chrs,unique(SNPpos[,1]))
#   
#   Tumor_LogR = Tumor_LogR[,c(-1,-2),drop=F]
#   Tumor_BAF = Tumor_BAF[,c(-1,-2),drop=F]
#   # make sure it is all converted to numerical values
#   for (cc in 1:dim(Tumor_LogR)[2]) {
#     Tumor_LogR[,cc]=as.numeric(as.vector(Tumor_LogR[,cc]))
#     Tumor_BAF[,cc]=as.numeric(as.vector(Tumor_BAF[,cc]))
#   }
#   if(!is.null(Germline_LogR_file)) {
#     Germline_LogR = Germline_LogR[,c(-1,-2),drop=F]
#     Germline_BAF = Germline_BAF[,c(-1,-2),drop=F]
#     for (cc in 1:dim(Germline_LogR)[2]) {
#       Germline_LogR[,cc]=as.numeric(as.vector(Germline_LogR[,cc]))
#       Germline_BAF[,cc]=as.numeric(as.vector(Germline_BAF[,cc]))
#     }
#   }
#   
#   # sort all data by genomic position
#   last = 0;
#   ch = list();
#   SNPorder = vector(length=dim(SNPpos)[1])
#   for (i in 1:length(chrs)) {
#     chrke = SNPpos[SNPpos[,1]==chrs[i],]
#     chrpos = chrke[,2]
#     names(chrpos) = rownames(chrke)
#     chrpos = sort(chrpos)
#     ch[[i]] = (last+1):(last+length(chrpos))  
#     SNPorder[ch[[i]]] = names(chrpos)
#     last = last+length(chrpos)
#   }
#   SNPpos = SNPpos[SNPorder,]
#   Tumor_LogR=Tumor_LogR[SNPorder,,drop=F]
#   Tumor_BAF=Tumor_BAF[SNPorder,,drop=F]
#   
#   if(!is.null(Germline_LogR_file)) {
#     Germline_LogR = Germline_LogR[SNPorder,,drop=F]
#     Germline_BAF = Germline_BAF[SNPorder,,drop=F]
#   }
#   
#   # split the genome into distinct parts to be used for segmentation (e.g. chromosome arms, parts of genome between gaps in array design)
#   print.noquote("Splitting genome in distinct chunks...")
#   chr = split_genome(SNPpos)
#   
#   if (is.null(gender)) {
#     gender = rep("XX",dim(Tumor_LogR)[2])
#   }
#   return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF, 
#               Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL, 
#               Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF, 
#               SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs, 
#               samples = colnames(Tumor_LogR), gender = gender, 
#               sexchromosomes = sexchromosomes,
#               failedarrays = NULL))
# }


#' Parse the reference info file
#' @param snp6_reference_info_file A SNP6 reference info master file
#' @noRd
parseSNP6refFile = function(snp6_reference_info_file) {
  return(read.table(snp6_reference_info_file, header=T, stringsAsFactors=F))
}

#' Transform cel files into BAF and LogR
#'
#' This function takes a cel file from a tumour and a matched normal and
#' extracts the BAF and LogR, which is saved into a single file. The \code{gc.correct}
#' function can read that file and transforms it into separate BAF and LogR files that
#' both Battenberg and ASCAT can use.
#' @param normal_cel_file String that points to the cel file containing the matched normal data
#' @param tumour_cel_file String that points to the cel file containing the tumour data
#' @param output_file String where the BAF and LogR should be written
#' @param snp6_reference_info_file String to the SNP6 reference info file that comes with Battenberg SNP6
#' @param apt.probeset.genotype.exe Path to the apt.probeset.genotype executable (Default $PATH)
#' @param apt.probeset.summarize.exe Path to the apt.probeset.summarize executable (Default $PATH)
#' @param norm.geno.clust.exe Path to the normalize_affy_geno_cluster.pl script (Default $PATH)
#' @author sd11
#' @export
cel2baf.logr = function(normal_cel_file, tumour_cel_file, output_file, snp6_reference_info_file, apt.probeset.genotype.exe="apt-probeset-genotype", apt.probeset.summarize.exe="apt-probeset-summarize", norm.geno.clust.exe="normalize_affy_geno_cluster.pl") {
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
  cmd = paste(norm.geno.clust.exe, UNM_NORMALS, "quant-norm.pm-only.med-polish.expr.summary.txt", "-locfile", LOCFILE, "-out", output_file)
  print(cmd)
  system(cmd, wait=T)
}

#' Correct the LogR estimates for GC content
#' 
#' This function performs GC correction of the LogR
#' data. Sometimes a wave pattern is observed there
#' that correlates with GC content. Internally it uses
#' the ASCAT gc correction function.
#' @param samplename Name of the sample to be used to name columns
#' @param infile.logr.baf String that points to the raw combined BAF and LogR file that is the result of \code{cel2baf.logr}
#' @param outfile.tumor.LogR The filename of the file where the tumour LogR will be written
#' @param outfile.tumor.BAF The filename of the file where the tumour BAF will be written
#' @param outfile.normal.LogR The filename of the file where the normal LogR will be written
#' @param outfile.normal.BAF The filename of the file where the normal BAF will be written
#' @param outfile.probeBAF The filename of the file where the probe ids and their BAF will be saved
#' @param snp6_reference_info_file String to the SNP6 reference info file that comes with Battenberg SNP6
#' @param chr_names A vector of chromosome names that are to be used
#' @param birdseed_report_file Name of the birdseed output file. This is a temp output file of one of the internally called functions of which the name cannot be defined. Don't change this parameter. (Default birdseed.report.txt)
#' @author sd11
#' @export
gc.correct = function(samplename, infile.logr.baf, outfile.tumor.LogR, outfile.tumor.BAF, outfile.normal.LogR, outfile.normal.BAF, outfile.probeBAF, snp6_reference_info_file, chr_names, birdseed_report_file="birdseed.report.txt") {
  # Read in needed reference files
  ref.files = parseSNP6refFile(snp6_reference_info_file)
  SNP_POS_REF = ref.files[ref.files$variable == "SNP_POS",]$reference_file
  GC_SNP6 = ref.files[ref.files$variable == "GC_SNP6",]$reference_file
  
  lrrbaf = read.table(infile.logr.baf, header=T, sep="\t", row.names=1, stringsAsFactors=F)
  SNPpos = read.table(SNP_POS_REF, header=T, sep="\t", row.names=1, stringsAsFactors=F)
  
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
  
  write.table(cbind(SNPpos,Tumor_BAF), paste(outfile.tumor.BAF, "_noGCcorr.txt", sep=""), sep="\t", row.names=T, quote=F)
  write.table(cbind(SNPpos,Normal_BAF), paste(outfile.normal.BAF, "_noGCcorr.txt", sep=""), sep="\t", row.names=T, quote=F)
  
  # read into ASCAT and make GC corrected input:
  write.table(cbind(SNPpos,Tumor_LogR), paste(outfile.tumor.LogR, "_noGCcorr.txt", sep=""), sep="\t", row.names=T, quote=F)
  write.table(cbind(SNPpos,Normal_LogR), paste(outfile.normal.LogR, "_noGCcorr.txt", sep=""), sep="\t", row.names=T, quote=F)
  
  # ======================================= above previous prepareGCcorrect, below runGCcorrect ==============================================
  
  # TODO: This must be a dapted to not hardcode the chromosome names
  gender <- read.table(birdseed_report_file, sep="\t", skip=66, header=T)
  sex <- as.vector(gender[,"computed_gender"])
  sex[sex == "female"] <- "XX"
  sex[sex == "male"] <- "XY"
  sex[sex == "unknown"] <- NA
  
  ascat.bc <- ASCAT::ascat.loadData(paste(outfile.tumor.LogR, "_noGCcorr.txt", sep=""), paste(outfile.tumor.BAF, "_noGCcorr.txt", sep=""),paste(outfile.normal.LogR, "_noGCcorr.txt", sep=""), paste(outfile.normal.BAF, "_noGCcorr.txt", sep=""), chrs=chr_names, gender=sex)
  ASCAT::ascat.plotRawData(ascat.bc)
  ascat.bc <- ASCAT::ascat.GCcorrect(ascat.bc, GC_SNP6)

  # Make sure the right column names are added here, because these are expected by fitcopynumber
  colnames(ascat.bc$SNPpos) = c("Chromosome", "Position")

  # Determine SNPs with BAF between 0.3-0.7 from normal => these are supposed to be heterozygous
  is.het = (ascat.bc$Germline_BAF >= 0.3 & ascat.bc$Germline_BAF <= 0.7)
  dat = cbind(ascat.bc$SNPpos, round(ascat.bc$Germline_LogR, 4))
  dat = dat[which(is.het),]
  colnames(dat) = c("Chromosome", "Position", samplename)
  write.table(dat, file=outfile.normal.LogR, row.names=F, quote=F, sep="\t")

  select = !is.na(ascat.bc$Germline_BAF)
  dat = cbind(ascat.bc$SNPpos, round(ascat.bc$Germline_BAF, 4))
  colnames(dat) = c("Chromosome", "Position", samplename)
  write.table(dat[which(select),], file=outfile.normal.BAF, row.names=F, quote=F, sep="\t")

  # Save the probe ids plus their BAF for only the germline heterozygous mutations
  select = !is.na(ascat.bc$Tumor_BAF)
  dat = cbind(row.names(ascat.bc$SNPpos), ascat.bc$Tumor_BAF)
  dat = dat[which(select & is.het),]
  write.table(dat, file=outfile.probeBAF, row.names=F, quote=F, col.names=F, sep="\t")

  # Save tumour BAF and LogR directly. Include homozygous SNPs here.
  dat = cbind(ascat.bc$SNPpos, round(ascat.bc$Tumor_BAF, 4))
  dat = dat[which(select),]
  colnames(dat) = c("Chromosome", "Position", samplename)
  write.table(dat, file=outfile.tumor.BAF, row.names=F, quote=F, sep="\t")

  select = !is.na(ascat.bc$Tumor_LogR)
  dat = cbind(ascat.bc$SNPpos, round(ascat.bc$Tumor_LogR, 4))
  dat = dat[which(select),]
  colnames(dat) = c("Chromosome", "Position", samplename)
  write.table(dat, file=outfile.tumor.LogR, row.names=F, quote=F, sep="\t")
}


#' Prepares data for impute
#'
#' The raw BAF and LogR data have been dumped into separate files. Now the data
#' needs to be prepared to go into Impute2, which is essentially morphing it into
#' the correct format. This function does that per chromosome and can therefore
#' be run in parallel for each chromosome.
#' @param infile.germlineBAF Germline BAF file generated by \code{cel2baf.logr}
#' @param infile.tumourBAF Tumour BAF file generated by \code{cel2baf.logr}
#' @param outFileStart Prefix of the filenames where the Impute2 input will be written. These will be extended with the chromosome
#' @param chrom Char with the chromosome for which an Impute2 file is produced
#' @param chr_names A vector of chromosome names that can be considered. This vector can just contain the chromosome for which the Impute2 file is produced, but can contain all chromosomes.
#' @param problemLociFile A string that points to a file with problematic loci that should be removed from the data
#' @param snp6_reference_info_file String to the SNP6 reference info file that comes with Battenberg SNP6
#' @param imputeinfofile String to the impute 1000 genomes reference info file that comes with Battenberg
#' @param is.male Boolean that is True if the donor is male, False when female
#' @param heterozygousFilter BAF cutoff for calling homozygous SNPs
#' @author dw9 jd
#' @export
generate.impute.input.snp6 = function(infile.germlineBAF, infile.tumourBAF, outFileStart, chrom, chr_names, problemLociFile, snp6_reference_info_file, imputeinfofile, is.male, heterozygousFilter="none") {
  # Obtain pointer to SNP6 specific reference file
  ref.files = parseSNP6refFile(snp6_reference_info_file)
  ANNO_FILE = ref.files[ref.files$variable == "ANNO_FILE",]$reference_file
  
  # Read in the 1000 genomes reference file paths for the specified chrom
  impute.info = parse.imputeinfofile(imputeinfofile, is.male, chrom=chrom)
  
  # Read in the known SNP locations from the 1000 genomes reference files
  known_SNPs = read.table(impute.info$impute_legend[1], sep=" ", header=T)
  if(nrow(impute.info)>1){
    for(r in 2:nrow(impute.info)){
      known_SNPs = rbind(known_SNPs, read.table(impute.info$impute_legend[r], sep=" ", header=T))
    }
  }
  
  outfile=paste(outFileStart,chrom,".txt",sep="")

  known_SNPs[,3]=factor(known_SNPs[,3],levels=c("A","C","G","T"))
  known_SNPs[,4]=factor(known_SNPs[,4],levels=c("A","C","G","T"))

  print(head(known_SNPs))
  print(dim(known_SNPs))
  chr_name = chrom
  
  # filter out bad SNPs (streaks in BAF)
  if((problemLociFile !="NA") & (!is.na(problemLociFile)))
  {
    problemSNPs=read.table(problemLociFile,header=T,sep="\t")
    problemSNPs=problemSNPs$Pos[problemSNPs$Chr==chr_name]
    badIndices=match(known_SNPs[,2],problemSNPs)
    known_SNPs = known_SNPs[is.na(badIndices),]
    print(paste("badIndices lengths=",length(badIndices),",",sum(is.na(badIndices)),sep=""))
  }
  
  knownSNP6data=read.csv(ANNO_FILE,comment.char="#",header=T,row.names=NULL,stringsAsFactors=F)
  knownSNP6data=knownSNP6data[knownSNP6data$Chromosome==chr_name,]
  print(paste("first column=",names(knownSNP6data)[1],sep=""))
  print(paste("first known datum=",knownSNP6data[1,1],sep=""))
  
  # adjust for strand
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
  
  # Read in the BAFs and see which 1000 genomes SNPs are covered
  germline_snp_data = read.table(infile.germlineBAF,sep="\t",header=T, stringsAsFactors=F) #[,3,drop=F]
  germline_snp_data = germline_snp_data[germline_snp_data[,1]==chr_name,]
  tumour_snp_data = read.table(infile.tumourBAF,sep="\t",header=T, stringsAsFactors=F) #[,3,drop=F]
  tumour_snp_data = tumour_snp_data[tumour_snp_data[,1]==chr_name,]
  # snp_matches = match(rownames(germline_snp_data), rownames(tumour_snp_data))
  snp_matches = match(germline_snp_data[,2], tumour_snp_data[,2])
  snp_data = na.omit(cbind(nBAF = germline_snp_data[,3], tBAF = tumour_snp_data[snp_matches,3]))

  print(paste("first datum=",rownames(snp_data[1,]),sep=""))

  #indices = match(rownames(snp_data),knownSNP6data$Probe.Set.ID)
  indices = match(germline_snp_data[,2], knownSNP6data$Physical.Position)
  if(sum(!is.na(indices))==0){
    print("Did not find any positional matches of the provided data to the reference")
    # indices = match(rownames(snp_data),knownSNP6data$dbSNP.RS.ID)
    q(save="no", status=1)
  }

  print(paste("found SNPs=",sum(!is.na(indices)),sep=""))
  print(paste("class=",class(knownSNP6data$Physical.Position),sep=""))
  matched.info = cbind(knownSNP6data[indices[!is.na(indices)],c("Physical.Position","Allele.A","Allele.B")],snp_data[!is.na(indices),1:2])
  print(paste("class2=",class(matched.info[,1]),sep=""))

  print(paste("first row of matched.info=",paste(matched.info[1,],sep=","),sep=""))
  print(paste("first Allele.A=",matched.info$Allele.A[1],sep=""))
  print(paste("first Allele.B=",matched.info$Allele.B[1],sep=""))
  
  print(paste("class 1a =",class(known_SNPs[,2]),sep=""))
  print(paste("class 2a =",class(as.numeric(known_SNPs[,2])),sep=""))

  indices2 = match(matched.info[,1],known_SNPs[,2])
  combined.info = na.omit(cbind(matched.info[!is.na(indices2),],known_SNPs[indices2[!is.na(indices2)],1:4]))
  print(paste("first row of combined.info=",paste(combined.info[1,],sep=","),sep=""))
  lev2 = levels(combined.info[,2])
  print(paste("levels[2]=",paste(lev2,sep=","),sep=""))
  lev3 = levels(combined.info[,3])
  print(paste("levels[3]=",paste(lev3,sep=","),sep=""))
  lev8 = levels(combined.info[,8])
  print(paste("levels[8]=",paste(lev8,sep=","),sep=""))
  lev9 = levels(combined.info[,9])
  print(paste("levels[9]=",paste(lev9,sep=","),sep=""))

  combined.info1 = combined.info[(combined.info[,2]==combined.info[,8] & combined.info[,3]==combined.info[,9]),]
  #alleles are reversed
  combined.info2 = cbind(combined.info[(combined.info[,2]==combined.info[,9] & combined.info[,3]==combined.info[,8]),1:3],1.0-combined.info[(combined.info[,2]==combined.info[,9] & combined.info[,3]==combined.info[,8]),4:5],combined.info[(combined.info[,2]==combined.info[,9] & combined.info[,3]==combined.info[,8]),6:9])
  names(combined.info2) = names(combined.info1)

  all.info = rbind(combined.info1,combined.info2)
  print(paste("norows all.info=",nrow(all.info),sep=""))

  all.info = all.info[order(as.numeric(all.info[,1])),]

  is.het = (all.info[,4] >= 0.3 & all.info[,4] <= 0.7)
  names(all.info)[5]="allele.frequency"
  write.csv(all.info[is.het,-4], file=paste(outFileStart,chrom,"_withAlleleFreq.csv",sep=""), quote=F, row.names=F)

  out.data = data.frame()
  if (heterozygousFilter!="none") {
    # Set the minimum level to use for calling homozygous SNPs
    minBaf = min(heterozygousFilter, 1.0-heterozygousFilter)
    maxBaf = max(heterozygousFilter, 1.0-heterozygousFilter)

    is.hom.ref = (all.info[,4] <= minBaf)
    is.hom.alt = (all.info[,4] >= maxBaf)

    # Obtain genotypes that impute2 is able to understand
    genotypes = array(0,c(nrow(all.info),3))
    genotypes[is.hom.ref,1] = 1
    genotypes[is.het,2] = 1
    genotypes[is.hom.alt,3] = 1
    is.genotyped = (is.het | is.hom.ref | is.hom.alt)

    snp.names = paste("snp",1:sum(is.genotyped),sep="")
    out.data = cbind(snp.names,all.info[is.genotyped,6:9], genotypes[is.genotyped,])
  } else {
    snp.names = paste("snp",1:sum(is.het),sep="")
    out.data = cbind(snp.names,all.info[is.het,6:9], matrix(data=c(0,1,0), nrow=sum(is.het), ncol=3, byrow=T))
  }
  write.table(out.data,file=outfile,row.names=F,col.names=F,quote=F)
  
  if (chrom=='chrX') {
    sample.g.file = paste(outFileStart,"sample_g.txt",sep="")
    sample_g_data = data.frame(ID_1=c(0,"INDIVI1"),ID_2=c(0,"INDIVI1"),missing=c(0,0),sex=c("D",2))
    write.table(sample_g_data, file=sample.g.file, row.names=F, col.names=T, quote=F)
  }
}

#' Infer the gender using the birdseed report file
#' @param birdseed_report_file The birdseed report file
#' @export
infer_gender_birdseed = function(birdseed_report_file) {
  z = read.table(birdseed_report_file, header=T)
  return(as.character(z$em.cluster.chrX.het.contrast_gender))
}


#' Prepare SNP6 data for haplotype construction
#' 
#' This function performs part of the Battenberg SNP6 pipeline: Extract BAF and logR from the CEL files 
#' and performing GC content correction.
#'
#' @param tumour_cel_file Full path to a CEL file containing the tumour raw data
#' @param normal_cel_file Full path to a CEL file containing the normal raw data
#' @param tumourname Identifier to be used for tumour output files
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param snp6_reference_info_file Full path to the SNP6 reference info file
#' @param apt.probeset.genotype.exe Full path to the apt.probeset.genotype executable (Default: expected in $PATH)
#' @param apt.probeset.summarize.exe Full path to the apt.probeset.summarize executable (Default: expected in $PATH)
#' @param norm.geno.clust.exe  Full path to the norm.geno.clust.exe executable (Default: expected in $PATH)
#' @param birdseed_report_file Name of the birdseed output file. This is a temp output file of one of the internally called functions of which the name cannot be defined. Don't change this parameter. (Default: birdseed.report.txt)
#' @author sd11
#' @export
prepare_snp6 = function(tumour_cel_file, normal_cel_file, tumourname, chrom_names, 
                        snp6_reference_info_file, apt.probeset.genotype.exe="apt-probeset-genotype", 
                        apt.probeset.summarize.exe="apt-probeset-summarize", norm.geno.clust.exe="normalize_affy_geno_cluster.pl",
                        birdseed_report_file="birdseed.report.txt") {
  
  # Extract the LogR and BAF from both tumour and normal cel files.
  cel2baf.logr(normal_cel_file=normal_cel_file,
               tumour_cel_file=tumour_cel_file,
               output_file=paste(tumourname, "_lrr_baf.txt", sep=""),
               snp6_reference_info_file=snp6_reference_info_file,
               apt.probeset.genotype.exe=apt.probeset.genotype.exe,
               apt.probeset.summarize.exe=apt.probeset.summarize.exe,
               norm.geno.clust.exe=norm.geno.clust.exe)

  gc.correct(samplename=tumourname,
             infile.logr.baf=paste(tumourname, "_lrr_baf.txt", sep=""),
             outfile.tumor.LogR=paste(tumourname, "_mutantLogR.tab", sep=""),
             outfile.tumor.BAF=paste(tumourname, "_mutantBAF.tab", sep=""),
             outfile.normal.LogR=paste(tumourname, "_germlineLogR.tab", sep=""),
             outfile.normal.BAF=paste(tumourname, "_germlineBAF.tab", sep=""),
             outfile.probeBAF=paste(tumourname, "_probeBAF.txt", sep=""),
             snp6_reference_info_file=snp6_reference_info_file,
             birdseed_report_file=birdseed_report_file,
             chr_names=chrom_names)
  
}
