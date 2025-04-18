% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmentation.R
\name{segment.baf.phased.multisample}
\alias{segment.baf.phased.multisample}
\title{Segment BAF, with the possible inclusion of structural variant breakpoints}
\usage{
segment.baf.phased.multisample(
  samplename,
  inputfile,
  outputfile,
  prior_breakpoints_file = NULL,
  gamma = 10,
  calc_seg_baf_option = 3,
  GENOMEBUILD
)
}
\arguments{
\item{samplename}{Name of the sample, which is used to name output figures}

\item{inputfile}{String that points to the output from the \code{combine.baf.files} function. This contains the phased SNPs with their BAF values}

\item{outputfile}{String where the segmentation output will be written}

\item{prior_breakpoints_file}{String that points to a file with prior breakpoints (from SVs for example) with chromosome and position columns (Default: NULL)}

\item{gamma}{The gamma parameter controls the size of the penalty of starting a new segment during segmentation. It is therefore the key parameter for controlling the number of segments (Default 10)}

\item{calc_seg_baf_option}{Various options to recalculate the BAF of a segment. Options are: 1 - median, 2 - mean, 3 - ifelse median==0 or 1, median, mean. (Default: 3)}

\item{GENOMEBUILD}{Genome build upon which the 1000G SNP coordinates were obtained}
}
\description{
This function breaks the genome up into chromosomes, possibly further when SV breakpoints
are provided, and runs PCF on each to segment the chromosomes independently.
}
\author{
jdemeul, sd11
}
