#!/bin/bash
samplename=$1

tar czf ${samplename}_battenberg_segmentation.tar.gz *egmented*txt *_segment_chr*png *RAFseg*png && rm *egmented*txt *_segment_chr*png *RAFseg*png
tar czf ${samplename}_battenberg_haplotyping.tar.gz *heterozygousMutBAFs_haplotyped.txt *_heterozygousData.png *impute_input* *_allHaplotypeInfo.txt && rm *heterozygousMutBAFs_haplotyped.txt *_heterozygousData.png *impute_input* *_allHaplotypeInfo.txt
tar czf ${samplename}_battenberg_alleleCounts.tar.gz *alleleFreq*chr* && rm *alleleFreq*chr*
if `ls | grep "mutant\|germline" | grep txt > /dev/null`; then
  # SNP6
  tar czf ${samplename}_battenberg_raw_baf_logr.tar.gz *mutant*txt *germline*txt *mutant*tab *germline*tab *lleleCounts.tab && rm *mutant*txt *germline*txt *mutant*tab *germline*tab *lleleCounts.tab
else 
  # WGS
  tar czf ${samplename}_battenberg_raw_baf_logr.tar.gz *mutant*tab *normal*tab *lleleCounts.tab && rm *mutant*tab *normal*tab *lleleCounts.tab
fi
tar czf ${samplename}_battenberg_other.tar.gz *solution_status.txt *GCwindowCorrelations* *segment_masking_details.txt sample_g.txt *subclones_1.txt *second_distance.png ${samplename}_nonroundedprofile.png ${samplename}_copynumberprofile.png && rm *solution_status.txt *GCwindowCorrelations* *segment_masking_details.txt sample_g.txt *subclones_1.txt *second_distance.png ${samplename}_nonroundedprofile.png ${samplename}_copynumberprofile.png
tar czf ${samplename}_battenberg_copynumber.tar.gz *subclones_chr*png *subclones.txt *_second_*png *rho_and_psi.txt *tumour.png *alleleratio.png *BattenbergProfile*png *coverage.png *distance.png *germline.png *totalcn_chrom_plot.png *cellularity_ploidy.txt *refit_suggestion.txt *Rout && rm *subclones_chr*png *subclones.txt *_second_*png *rho_and_psi.txt *tumour.png *alleleratio.png *BattenbergProfile*png *coverage.png *distance.png *germline.png *totalcn_chrom_plot.png *cellularity_ploidy.txt *refit_suggestion.txt *Rout