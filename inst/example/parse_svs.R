library(VariantAnnotation)

#' Parses BRASS SV calls into a dataframe with a line for each SV and two columns: chromosome and position
parse_brass_svs = function(vcffile, outfile, ref_genome="hg19") {
  svs = parse_svs_1(vcffile, ref_genome=ref_genome)
  write_svs(svs, outfile)
  return(svs)
}

#' Parses ICGC consensus SV calls into a dataframe with a line for each SV and two columns: chromosome and position
parse_icgc_consensus_svs = function(vcffile, outfile, ref_genome="hg19") {
  svs = parse_svs_1(vcffile, ref_genome=ref_genome)
  write_svs(svs, outfile)
  return(svs)
}

#' Helper function that writes the given SVs to file
write_svs = function(svs, filename) {
  write.table(svs, file=filename, quote=F, row.names=F, sep="\t")
}

#' Helper function that works on cases where SVs have been encoded as such:
#' 
#' #CHROM  POS     ID      REF
#' 1       123        1     X       X]6:578]   
#' 1       234        2     X       ]1:280]YYX 
#' 1       280        3     Z       ZYY[1:234[
parse_svs_1 = function(vcffile, ref_genome="hg19") {
  svs = readVcf(vcffile, genome=ref_genome)
  output = data.frame(chromosome=seqnames(svs), position=start(svs))
  endpoints = alt(svs)
  endpoints = lapply(endpoints, function(x) {
    if (grepl("[", x, fixed=T)) {
      chrompos = unlist(strsplit(x, "[", fixed=T))[2]
    } else if (grepl("]", x, fixed=T)) {
      chrompos = unlist(strsplit(x, "]", fixed=T))[2]
    } else {
      chrompos = NA
    }
    chrompos_split = unlist(strsplit(chrompos, ":", fixed=T))
    return(data.frame(chromosome=chrompos_split[1], position=as.numeric(chrompos_split[2])))
  })
  endpoints = do.call(rbind, endpoints)
  output = rbind(output, endpoints)
  output = with(output, output[order(chromosome, position),])
  output$chromosome = as.character(output$chromosome)
  output = unique(output)
  output = output[with(output, order(chromosome, position)),]
  return(output)
}

#' Helper function that works on cases where SVs have been encoded as such:
#' 
#' #CHROM  POS     ...      INFO
#' 1       123     ...     ...CHR2=1;END=143274758...
#' 1       234     ...     ...CHR2=1;END=143274758...
#' 1       280     ...     ...CHR2=1;END=143274758...
parse_svs_2 = function(vcffile, ref_genome="hg19") {
  v = readVcf(vcffile, ref_genome)
  output = data.frame(chromosome=seqnames(v), position=start(v))
  output = rbind(output, data.frame(chromosome=info(v)$CHR2, position=info(v)$END))
  return(output)
}

parse_delly_svs = function(vcffile, outfile, ref_genome="hg19") {
  svs = parse_svs_2(vcffile, ref_genome)
  write_svs(svs, outfile)
  return(svs)
}





