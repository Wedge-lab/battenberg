# helper function to split the genome into parts
split_genome = function(SNPpos) {

  # look for gaps of more than 1Mb and chromosome borders
  holesOver1Mb = which(diff(SNPpos[,2])>=1000000)+1
  chrBorders = which(diff(as.numeric(SNPpos[,1]))!=0)+1
  holes = unique(sort(c(holesOver1Mb,chrBorders)))

  # find which segments are too small
  joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)

  # if it's the first or last segment, just join to the one next to it, irrespective of chromosome and positions
  while (1 %in% joincandidates) {
    holes=holes[-1]
    joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  }
  while ((length(holes)+1) %in% joincandidates) {
    holes=holes[-length(holes)]
    joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  }
 
  while(length(joincandidates)!=0) {
    # the while loop is because after joining, segments may still be too small..

    startseg = c(1,holes)
    endseg = c(holes-1,dim(SNPpos)[1])

    # for each segment that is too short, see if it has the same chromosome as the segments before and after
    # the next always works because neither the first or the last segment is in joincandidates now
    previoussamechr = SNPpos[endseg[joincandidates-1],1]==SNPpos[startseg[joincandidates],1] 
    nextsamechr = SNPpos[endseg[joincandidates],1]==SNPpos[startseg[joincandidates+1],1]

    distanceprevious = SNPpos[startseg[joincandidates],2]-SNPpos[endseg[joincandidates-1],2]
    distancenext = SNPpos[startseg[joincandidates+1],2]-SNPpos[endseg[joincandidates],2]

    # if both the same, decide based on distance, otherwise if one the same, take the other, if none, just take one.
    joins = ifelse(previoussamechr&nextsamechr, 
                   ifelse(distanceprevious>distancenext, joincandidates, joincandidates-1),
                   ifelse(nextsamechr, joincandidates, joincandidates-1))

    holes=holes[-joins]

    joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  }
  # if two neighboring segments are selected, this may make bigger segments then absolutely necessary, but I'm sure this is no problem.

  startseg = c(1,holes)
  endseg = c(holes-1,dim(SNPpos)[1])

  chr=list()
  for (i in 1:length(startseg)) {
    chr[[i]]=startseg[i]:endseg[i]
  }
  
  return(chr)
}

# function to read in SNP array data
# input: filenames of tumor LogR, tumor BAF, germline LogR and germline BAF data
# germline data files can be NULL - in that case these are not read in
# output: a data structure containing:
# 1. Tumor_LogR data matrix
# 2. Tumor_BAF data matrix
# 3. Tumor_LogR_segmented: placeholder, NULL
# 4. Tumor_BAF_segmented: placeholder, NULL
# 5. Germline_LogR data matrix
# 6. Germline_BAF data matrix
# 7. SNPpos: position of all SNPs
# 8. ch: a list containing vectors with the indices for each chromosome (e.g. Tumor_LogR[ch[[13]],] will output the Tumor_LogR data of chromosome 13
# 9. chr: a list containing vectors with the indices for each distinct part that can be segmented separately (e.g. chromosome arm, stretch of DNA between gaps in the array design)
ascat.loadData = function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, Germline_BAF_file = NULL, chrs = c(1:22,"X"), Tumor_counts_file = NULL, Germline_counts_file = NULL) {

  Tumor_counts = NULL
  Germline_counts = NULL

  # read in SNP array data files
  print.noquote("Reading Tumor LogR data...")
  Tumor_LogR <- read.table(Tumor_LogR_file, header=T, row.names=1, comment.char="", sep = "\t")
  print.noquote("Reading Tumor BAF data...")
  Tumor_BAF <- read.table(Tumor_BAF_file, header=T, row.names=1, comment.char="", sep = "\t")
  if(!is.null(Tumor_counts_file)){
	  Tumor_counts <- read.table(Tumor_counts_file, header=T, row.names=1, comment.char="", sep = "\t")
  }

  Germline_LogR = NULL
  Germline_BAF = NULL
  if(!is.null(Germline_LogR_file)) {
    print.noquote("Reading Germline LogR data...")
    Germline_LogR <- read.table(Germline_LogR_file, header=T, row.names=1, comment.char="", sep = "\t")
    print.noquote("Reading Germline BAF data...")
    Germline_BAF <- read.table(Germline_BAF_file, header=T, row.names=1, comment.char="", sep = "\t")
	if(!is.null(Germline_counts_file)){
		  Germline_counts <- read.table(Germline_counts_file, header=T, row.names=1, comment.char="", sep = "\t")
	}    
  }

  # make SNPpos vector that contains genomic position for all SNPs and remove all data not on chromosome 1-22,X
  print.noquote("Registering SNP locations...")
  SNPpos <- Tumor_LogR[,1:2]
  SNPpos = SNPpos[SNPpos[,1]%in%chrs,]
  
	#print.noquote("OK1")

  Tumor_LogR = Tumor_LogR[rownames(SNPpos),c(-1,-2),drop=F]
  Tumor_BAF = Tumor_BAF[rownames(SNPpos),c(-1,-2),drop=F]
  if(!is.null(Tumor_counts_file)){
	Tumor_counts = Tumor_counts[rownames(SNPpos),]
  }
  
	#print.noquote("OK2")
	
  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[rownames(SNPpos),c(-1,-2),drop=F]
    Germline_BAF = Germline_BAF[rownames(SNPpos),c(-1,-2),drop=F]
    if(!is.null(Germline_counts_file)){
		Germline_counts = Germline_counts[rownames(SNPpos),]
    }    
  }

	#print.noquote("OK3")
	
  # sort all data by genomic position
  last = 0;
  ch = list();
  SNPorder = vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
	#print.noquote(paste("chr=",chrs[i],sep=""))
	  	
    chrke = SNPpos[SNPpos[,1]==chrs[i],]
    chrpos = chrke[,2]
    names(chrpos) = rownames(chrke)
    #print.noquote(paste("Sorting chromosome ",chrs[i],"...",sep=""))
    #time1=Sys.time()
    chrpos = sort(chrpos)
    #DCW if-else added 020611
    if(length(chrpos)>0){
		ch[[i]] = (last+1):(last+length(chrpos))  
		SNPorder[ch[[i]]] = names(chrpos)
		last = last+length(chrpos)
    }
    else{
    	ch[[i]] = array(0,0)
    }
    #time2=Sys.time()
    #print.noquote(paste("Finished sorting chromosome ",chrs[i],". Time taken= ",difftime(time2,time1),sep=""))
  }
  SNPpos = SNPpos[SNPorder,]
  Tumor_LogR=Tumor_LogR[SNPorder,,drop=F]
  Tumor_BAF=Tumor_BAF[SNPorder,,drop=F]
  if(!is.null(Tumor_counts_file)){
	  Tumor_counts=Tumor_counts[SNPorder,]
  }
  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[SNPorder,,drop=F]
    Germline_BAF = Germline_BAF[SNPorder,,drop=F]
  	if(!is.null(Germline_counts_file)){
	  	Germline_counts=Germline_counts[SNPorder,]
  	}    
  }

  # write final SNP positions to file
  #write.table(SNPpos,"SNPpos.txt", sep="\t")

  # split the genome into distinct parts to be used for segmentation (e.g. chromosome arms, parts of genome between gaps in array design)
  print.noquote("Splitting genome in distinct chunks...")
  chrom = split_genome(SNPpos)

  return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF, 
              Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL, 
              Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF, 
              Tumor_counts=Tumor_counts, Germline_counts=Germline_counts,
              SNPpos = SNPpos, ch = ch, chrom = chrom, chrs = chrs, 
              samples = colnames(Tumor_LogR)))
}




# plots SNP array data
# input: an ASCAT object (e.g. from ASCAT.loaddata and plots the SNP array data
# tumorfiles: start of filename for tumor data plots (no plotting if NULL)
# germlinefiles: start of filename for germline data plots (no plotting if NULL)
ascat.plotRawData = function(ASCATobj, tumorfiles = "Tumor", germlinefiles = "Germline", parentDir="./") {
  attach(ASCATobj)
  if(!is.null(tumorfiles)) {
    print.noquote("Plotting tumor data")
    for (i in 1:dim(Tumor_LogR)[2]) {
      png(filename = paste(parentDir,tumorfiles,samples[i],".png",sep=""), width = 5000, height = 2500, res = 500)
      par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(Tumor_LogR)[1]>20000,".",20))
      #par(mar = c(0.2,2,2,0.2), mfrow = c(2,1), cex = 1, cex.main=1.2, cex.axis = 0.8, pch = ifelse(dim(Tumor_LogR)[1]>50000,".",20))
      #par(mar = c(0.2,2,2,0.2), mfrow = c(2,1), cex = 1, cex.main=1.2, cex.axis = 0.8, pch = ".")
      plot(c(1,dim(Tumor_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(samples[i], ", tumor data, LogR", sep = ""), xlab = "", ylab = "")
      points(Tumor_LogR[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 0.8)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      plot(c(1,dim(Tumor_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(samples[i], ", tumor data, BAF", sep = ""), xlab = "", ylab = "")
      points(Tumor_BAF[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 0.8)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      dev.off()
    }
  }
  if(!is.null(germlinefiles) && !is.null(Germline_LogR)) {
    print.noquote("Plotting germline data")
    for (i in 1:dim(Germline_LogR)[2]) {
      png(filename = paste(parentDir,germlinefiles,samples[i],".png",sep=""), width = 5000, height = 2500, res = 500)
      par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(Tumor_LogR)[1]>20000,".",20))
      #par(mar = c(0.2,2,2,0.2), mfrow = c(2,1), cex = 1, cex.main=1.2, cex.axis = 0.8, pch = ifelse(dim(Tumor_LogR)[1]>50000,".",20))
      plot(c(1,dim(Germline_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(samples[i], ", germline data, LogR", sep = ""), xlab = "", ylab = "")
      points(Germline_LogR[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 0.8)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      plot(c(1,dim(Germline_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(samples[i], ", germline data, BAF", sep = ""), xlab = "", ylab = "")
      points(Germline_BAF[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ch)) {
        chrk = ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,chrs[j], pos = 1, cex = 0.8)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      dev.off()
    }
  }
  detach(ASCATobj)
}

# helping function to read segments:
make_seg_lr = function(r) {
  pcf_segments = numeric(0);
  index = 0;
  previousr = 1E10;
  for (i in 1:length(r)) {
    if (r[i] != previousr) {
      index=index+1;
      count=1;
    }
    else {
      count = count + 1;
    }
    pcf_segments[index] = count;
    previousr = r[i];
  }
  return(pcf_segments);
}

# function to make segments of constant LRR and BAF 
# this function is more general and does not depend on specifically ASPCF output
# it can also handle segmention performed on LRR and BAF separately
make_segments = function(r,b) {
  m = matrix(ncol = 2, nrow = length(b))
  m[,1] = r
  m[,2] = b
  m = as.matrix(na.omit(m))
  pcf_segments = matrix(ncol = 3, nrow = dim(m)[1])
  colnames(pcf_segments) = c("r","b","length");
  index = 0;
  previousb = -1;
  previousr = 1E10;
  for (i in 1:dim(m)[1]) {
    if (m[i,2] != previousb || m[i,1] != previousr) {
      index=index+1;
      count=1;
      pcf_segments[index, "r"] = m[i,1];
      pcf_segments[index, "b"] = m[i,2];
    }
    else {
      count = count + 1;
    }
    pcf_segments[index, "length"] = count;
    previousb = m[i,2];
    previousr = m[i,1];
  }
  pcf_segments = as.matrix(na.omit(pcf_segments))[,]
  return(pcf_segments);
}



# function to create the distance matrix (distance for a range of ploidy and tumor percentage values)
# input: segmented LRR and BAF and the value for gamma
#create_distance_matrix = function(segments, gamma) {
#DCW 180313 - user defined rho range allows forcing of rho to appropriate values (useful for very low rho)
create_distance_matrix = function(segments, gamma, possible.rho=seq(0.1,1.05,0.01)) {	
  s = segments
  psi_pos = seq(1,5.4,0.05) 
  rho_pos = possible.rho
  d = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) = psi_pos
  colnames(d) = rho_pos
  dmin = 1E20;
  for(i in 1:length(psi_pos)) {
    psi = psi_pos[i]
    for(j in 1:length(rho_pos)) {
      rho = rho_pos[j]
      nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
      nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
      # choose the minor allele
      nMinor = NULL
      if (sum(nA,na.rm=T) < sum(nB,na.rm=T)) {
        nMinor = nA
      }
      else {
        nMinor = nB
      }
      #d[i,j] = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1), na.rm=T)
      #DCW 180711 - try weighting BAF=0.5 equally with other points
      d[i,j] = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"], na.rm=T)
    }
  }
  return(d)
}



# the ASCAT main function
# lrr: (unsegmented) log R, in genomic sequence (all probes), with probe IDs
# baf: (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
# lrrsegmented: log R, segmented, in genomic sequence (all probes), with probe IDs
# bafsegmented: B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
# gamma: technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 % aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
# chromosomes: a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
# distancepng: if NA: distance is plotted, if filename is given, the plot is written to a .png file
# copynumberprofilespng: if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file
# nonroundedprofilepng: if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file
#runASCAT = function(lrr, baf, lrrsegmented, bafsegmented, chromosomes, distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA, gamma = 0.55, allow100percent,reliabilityFile=NA) {
runASCAT = function(lrr, baf, lrrsegmented, bafsegmented, chromosomes, distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA, gamma = 0.55, allow100percent, textOutput) {
	ch = chromosomes
  b = bafsegmented
  #r = lrrsegmented[names(bafsegmented)]
  #220213
  r = lrrsegmented

	
  library(RColorBrewer)

  s = make_segments(r,b)

  d = create_distance_matrix(s, gamma)

  if (!is.na(distancepng)) {
    png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
  }

  par(mar = c(5,5,0.5,0.5), cex=0.75, cex.lab=2, cex.axis=2)

  hmcol = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  image(log(d), col = hmcol, axes = F, xlab = "Ploidy", ylab = "Aberrant cell fraction")

  axis(1, at = seq(0, 4/4.4, by = 1/4.4), label = seq(1, 5, by = 1))
  axis(2, at = seq(0, 1/1.05, by = 1/3/1.05), label = seq(0.1, 1, by = 3/10))

  TheoretMaxdist = sum(rep(0.25,dim(s)[1]) * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1),na.rm=T)
  #DCW 180711 - try weighting BAF=0.5 equally with other points
  #TheoretMaxdist = sum(rep(0.25,dim(s)[1]) * s[,"length"],na.rm=T)

  nropt = 0
  localmin = NULL
  optima = list()
  for (i in 4:(dim(d)[1]-3)) {
    for (j in 4:(dim(d)[2]-3)) {
      m = d[i,j]
      seld = d[(i-3):(i+3),(j-3):(j+3)]
      seld[4,4] = max(seld)
      if(min(seld) > m) {
        psi = as.numeric(rownames(d)[i])
        rho = as.numeric(colnames(d)[j])
        nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
        nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
        
        # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
        ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
      
        percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
        perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
        # the next can happen if BAF is a flat line at 0.5 
        if (is.na(perczeroAbb)) {
          perczeroAbb = 0
        }

        goodnessOfFit = (1-m/TheoretMaxdist) * 100
		print(paste("ploidy=",ploidy,",rho=",rho,",goodness=",goodnessOfFit,",percentzero=",percentzero,", perczerAbb=",perczeroAbb,sep=""))
        if (ploidy > 1.6 & ploidy < 4.8 & rho >= 0.2 & goodnessOfFit > 80 & (percentzero > 0.01 | perczeroAbb > 0.1)) {
          nropt = nropt + 1
          optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
          localmin[nropt] = m
        }
      }     
    }
  }
  # if solutions with 100 % aberrant cell fraction should be allowed:
  # if there are no solutions, drop the conditions on regions with copy number zero, and include the borders (rho = 1) as well
  # this way, if there is another solution, this is still preferred, but these solutions aren't standardly eliminated
  if (allow100percent & nropt == 0) {
    #first, include borders
    cold = which(as.numeric(colnames(d))>1)
    d[,cold]=1E20
    for (i in 4:(dim(d)[1]-3)) {
      for (j in 4:(dim(d)[2]-3)) {
        m = d[i,j]
        seld = d[(i-3):(i+3),(j-3):(j+3)]
        seld[4,4] = max(seld)
        if(min(seld) > m) {
          psi = as.numeric(rownames(d)[i])
          rho = as.numeric(colnames(d)[j])
          nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
          nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
        
          # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
          ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
      
          percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
          perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
          # the next can happen if BAF is a flat line at 0.5 
          if (is.na(perczeroAbb)) {
            perczeroAbb = 0
          }

          goodnessOfFit = (1-m/TheoretMaxdist) * 100

		print(paste("ploidy etc",ploidy,rho,goodnessOfFit,sep=":"))

          #if (ploidy > 1.6 & ploidy < 4.8 & rho >= 0.2 & goodnessOfFit > 80) {
		  #180313 - don't check for minimum rho if possible.rho.values has been provided
		  if (is.na(possible.rho.values) & ploidy > 1.6 & ploidy < 4.8 & rho >= 0.2 & goodnessOfFit > 80 | !is.na(possible.rho.values) & ploidy > 1.6 & ploidy < 4.8 & goodnessOfFit > 65) {
            nropt = nropt + 1
            optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
            localmin[nropt] = m
          }
        }
      }
    }
  }

  if (nropt>0) {
    optlim = sort(localmin)[1]
    for (i in 1:length(optima)) {
      if(optima[[i]][1] == optlim) {  
        psi_opt1 = as.numeric(rownames(d)[optima[[i]][2]])
        rho_opt1 = as.numeric(colnames(d)[optima[[i]][3]])
        if(rho_opt1 > 1) {
          rho_opt1 = 1
        }
        ploidy_opt1 = optima[[i]][4]
        goodnessOfFit_opt1 = optima[[i]][5]
        points((psi_opt1-1)/4.4,(rho_opt1-0.1)/0.95,col="green",pch="X", cex = 2)
      }
    }
  }

  if (!is.na(distancepng)) {
    dev.off()
  }


  if(nropt>0) {

    if (!is.na(copynumberprofilespng)) {
      png(filename = copynumberprofilespng, width = 2000, height = 1000, res = 200)
    } 
    else {      
      windows(10,5)
    }

    par(mar = c(0.5,5,5,0.5), mfrow=c(2,1), cex = 0.4, cex.main=3, cex.axis = 2.5)

    rho = rho_opt1
    psi = psi_opt1
    
    write.table(data.frame(rho=rho,psi=psi),textOutput,sep="\t",row.names=F,quote=F)
    
    nAfull = (rho-1-(b-1)*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nBfull = (rho-1+b*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nA = pmax(round(nAfull),0)
    nB = pmax(round(nBfull),0)
    maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy_opt1),", aberrant cell fraction: ",sprintf("%2.0f",rho_opt1*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit_opt1),"%",sep="")
    plot(c(1,length(nAfull)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(nA-0.1,col="red",pch = "|")
    points(nB+0.1,col="green",pch = "|")
# don't ask me why, but the "|" sign is not centered, so the lines may need to be shifted..
    abline(v=0,lty=1,col="lightgrey")
    chrk_tot_len = 0
    
    print(names(lrr)[1:5])
    print(names(bafsegmented)[1:5])
    
    for (i in 1:length(ch)) {
      chrk = ch[[i]];
      print(chrk[1:5])
      chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
      print(chrk_hetero[1:5])
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      print(paste("drawing line for chr ",i," at ", tpos,sep=""))
      text(tpos,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)
      abline(v=vpos,lty=1,col="lightgrey")
    }
    
    rBacktransform = gamma*log((rho*(nA+nB)+(1-rho)*2)/((1-rho)*2+rho*psi),2)
    bBacktransform = (1-rho+rho*nB)/(2-2*rho+rho*(nA+nB))
    rConf = ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
    bConf = ifelse(bBacktransform!=0.5,pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))),NA)
 
    confidence = ifelse(is.na(rConf),bConf,ifelse(is.na(bConf),rConf,(rConf+bConf)/2))
    maintitle = paste("Aberration reliability score (%), average: ", sprintf("%2.0f",mean(confidence,na.rm=T)),"%",sep="")
    plot(c(1,length(nAfull)), c(0,100), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(confidence,col="blue",pch = "|")
    abline(v=0,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (i in 1:length(ch)) {
      chrk = ch[[i]];
      chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      print(paste("drawing line for chr ",i," at ", tpos,sep=""))
      text(tpos,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)    
      abline(v=vpos,lty=1,col="lightgrey")
    }
    

    if (!is.na(copynumberprofilespng)) {
      dev.off()
    }


    if (!is.na(nonroundedprofilepng)) {
      png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200)
    } 
    else {      
      windows(10,5)
    }

    par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)

    rho = rho_opt1
    psi = psi_opt1
    
    nAfull = (rho-1-(b-1)*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nBfull = (rho-1+b*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nA = pmax(round(nAfull),0)
    nB = pmax(round(nBfull),0)
    maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy_opt1),", aberrant cell fraction: ",sprintf("%2.0f",rho_opt1*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit_opt1),"%",sep="")
    plot(c(1,length(nAfull)), c(0,5), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
    points(nBfull,col="blue",pch = "|")
    points(nAfull+nBfull,col="purple",pch = "|")
# don't ask me why, but the "|" sign is not centered, so the lines may need to be shifted..
    abline(v=0,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (i in 1:length(ch)) {
      chrk = ch[[i]];
      chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk_hetero)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,5,ifelse(i<23,sprintf("%d",i),"X"), pos = 1, cex = 2)
      abline(v=vpos,lty=1,col="lightgrey")
    }
    
    if (!is.na(nonroundedprofilepng)) {
      dev.off()
    }
  }
}


