
library(optparse)
library(gtools)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input VCF file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output VCF file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

infile = opt$input
outfile = opt$output
genome = opt$genome

brass = read.table(infile, header=F, comment.char="#", stringsAsFactor=F)

# fetch  TRDS entry
trds_data = lapply(brass$V8, function(x) { r=unlist(strsplit(as.character(x), ";")); sel=grepl("TRDS", r); if (any(sel)) { unlist(strsplit(gsub("TRDS=", "", r[sel]), ",")) } else { NULL } })
brass$tumour_support = unlist(lapply(trds_data, length))

brass_filter = brass[brass$tumour_support > 0,]

# fetch second breakpoint, sometimes it is not mentioned as a first breakpoint
#second_chrpos = unlist(as.list(brass)$ALT))
second_chrpos = brass_filter$V5
second_chrpos = unlist(lapply(second_chrpos, function(x) {
       if (grepl("\\[", x)) { unlist(strsplit(x, "\\["))[2]
        } else { unlist(strsplit(x, "\\]"))[2] }}
       ))

brass_breakpoints = data.frame(chromosome=brass_filter$V1, position=brass_filter$V2, stringsAsFactors=F)
brass_breakpoints = rbind(brass_breakpoints, data.frame(chromosome=unlist(lapply(second_chrpos, function(x) unlist(strsplit(x, ":"))[1])), position=as.numeric(unlist(lapply(second_chrpos, function(x) unlist(strsplit(x, ":"))[2]))), stringsAsFactors=F))
brass_breakpoints = unique(brass_breakpoints)

# sort
brass_breakpoints_ordered = df = data.frame(matrix(ncol = 2, nrow = 0))
colnames(brass_breakpoints_ordered) = c("chromosome", "position")
for (chrom in gtools::mixedsort(unique(brass_breakpoints$chromosome))) {
	b_chrom = brass_breakpoints[brass_breakpoints$chromosome==chrom,]
	b_chrom = b_chrom[order(b_chrom$position),]
	brass_breakpoints_ordered = rbind(brass_breakpoints_ordered, b_chrom)
}
write.table(brass_breakpoints_ordered, file=outfile, quote=F, row.names=F, sep="\t")
