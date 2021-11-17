
####################################################################################################

#' A helper function to split the genome into parts
#' @param SNPpos A data.frame with a row for each SNP. First column is chromosome, second column position
#' @noRd
split_genome = function(SNPpos) {
  # look for gaps of more than 1Mb and chromosome borders
  holesOver1Mb = which(diff(SNPpos[,2])>=1000000)+1
  chrBorders = which(diff(as.numeric(factor(SNPpos[,1],levels=unique(SNPpos[,1]))))!=0)+1
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

####################################################################################################
#' Helper function that calculates a t-statistic
#' @noRd		
studentise <-function( sample_size, sample_mean, sample_SD, mu_pop )  # kjd 18-12-2013
{
	tvar = ( sample_mean - mu_pop ) * sqrt( sample_size ) / sample_SD
	
	return( tvar )
	
}

####################################################################################################
#' This function calculates a P-value, for a test where the null hypothesis is that
#' the sample was drawn from a Gaussian population with the specified mean "mu_pop".
#' @noRd
calc_Pvalue_t_twotailed <-function( sample_size, sample_mean, sample_SD, mu_pop, max_dist)  # kjd 18-12-2013
{
	tvar = ( sample_mean - mu_pop ) * sqrt( sample_size ) / sample_SD
	
	if( tvar < 0 )
	{
		lower_tail_prob = pt( tvar , df = sample_size - 1 , lower.tail = TRUE )
		
	}else
	{
		lower_tail_prob = 1 - pt( tvar , df = sample_size - 1 , lower.tail = TRUE )
		
	}
	
	pval = 2 * lower_tail_prob
	
	#DCW 250314
	if(abs(sample_mean - mu_pop)<max_dist){
		pval = 1
	}
	
	return( pval )
	
}

####################################################################################################
#' Helper function that calculates a binomial probability
#' @noRd
calc_binomial_prob <-function( sample_proportion, sample_size, pop_proportion ) # kjd 10-2-2014
{
	sample_count = round( sample_proportion * sample_size , 0 )
	
	if( sample_count < 0 ){ 
		sample_count = 0 
	}
	
	if( sample_count > sample_size ){ 
		sample_count = sample_size 
	}
	
	if( pop_proportion < 0 ){ 
		pop_proportion = 0 
	}
	
	if( pop_proportion > 1 ){ 
		pop_proportion = 1 
	}
	
	prob = dbinom( sample_count, sample_size, pop_proportion )
	
	return( prob )
	
}

####################################################################################################
#' This function calculates a log likelihood ratio where the two hypotheses are that
#' the tumour genome segment in question is "clonal".
#' The first hypothesis is the "best fit" model we can find. 
#' The second hypothesis is the "second best fit" model we can find. 
#' @noRd
calc_ln_likelihood_ratio <-function( LogR, BAFreq, BAF.length, BAF.size, BAF.mean, read_depth, rho, psi, gamma_param, maxdist_BAF ) # kjd 18-12-2013
{	
  pooled_BAF.size = read_depth * BAF.size
  
  # if we don't have a value for LogR, fill in 0
  if (is.na(LogR)) {
    LogR = 0
  }
  nMajor = (rho-1+BAFreq*psi*2^(LogR/gamma_param))/rho
  nMinor = (rho-1+(1-BAFreq)*psi*2^(LogR/gamma_param))/rho
  
  # to make sure we're always in a positive square:
  #if(nMajor < 0) {
  #	nMajor = 0.01
  #}
  #
  #if(nMinor < 0) {
  #	nMinor = 0.01
  #}
  #DCW - increase nMajor and nMinor together, to avoid impossible combinations (with negative subclonal fractions)
  if(nMinor<0 | is.na(nMinor)){
    if(BAFreq==1){
      #avoid calling infinite copy number
      nMajor = 1000
    }else{
      nMajor = nMajor + BAFreq * (0.01 - nMinor) / (1-BAFreq)
      if (nMajor<0) nMajor=1000
    }
    nMinor = 0.01		
  }

  if (!is.finite(nMajor)) {
    nMajor = 0.01
  }

  # Check if there is a viable solution
  if (!is.na(BAFreq)) { 
    nearest_edge = GetNearestCorners_bestOption( rho, psi, BAFreq, nMajor, nMinor ) # kjd 14-2-2014
    nMaj = nearest_edge$nMaj # kjd 14-2-2014
    nMin = nearest_edge$nMin # kjd 14-2-2014
  
  
    BAF_levels = (1-rho+rho*nMaj)/(2-2*rho+rho*(nMaj+nMin))
  
    index_vect = which( is.finite(BAF_levels) ) # kjd 14-2-2014
    BAF_levels = BAF_levels[ index_vect ] # kjd 14-2-2014
  
    if( length( BAF_levels ) > 1 ) # kjd 14-2-2014
    {
      likelihood_vect = sapply( BAF_levels , function(x){ calc_binomial_prob( BAF.mean, pooled_BAF.size, x ) } )
      likelihood_vect = sort( likelihood_vect, decreasing = TRUE )
    
      if( ( likelihood_vect[1] > 0 ) && ( likelihood_vect[2] > 0 ) )
      {
        ln_lratio = log( likelihood_vect[1] ) - log( likelihood_vect[2] )
      
      }else
      {
        ln_lratio = 0
      }
    
    }else
    {
      ln_lratio = 0
    }
  } else {
    ln_lratio = 0
  }
  
  return( ln_lratio )
  
}

####################################################################################################
		
#' Calculate a two tailed binomial p-value
#' @noRd
calc_Pvalue_binomial_twotailed <-function( sample_count, sample_size, pop_proportion ) # kjd 27-2-2014
{
	lower_tail_prob = pbinom( sample_count, sample_size, pop_proportion , lower.tail = TRUE )
	
	if( lower_tail_prob < 0.5 )
	{
		pval = 2 * lower_tail_prob
		
	}else
	{
		pval = 2 * ( 1 - lower_tail_prob )
		
	}
	
	return( pval )
	
}

####################################################################################################
#' Helper function that calculates a p-value for a set of BAF values summarised by their mean
#' TODO: this function is not used in Battenberg
#' @noRd
calc_BAF_Pvalue <-function( BAF.mean, pooled_BAF.size, maxdist_BAF, BAF_level ) # kjd 27-2-2014
{
	
	if( is.finite( BAF_level ) && pooled_BAF.size > 0 )
	{
		sample_size = round( pooled_BAF.size , 0 )
		sample_count = round( BAF.mean * pooled_BAF.size , 0 )
		
		if( sample_count < 0 ){ 
			sample_count = 0 
		}
		
		if( sample_count > sample_size ){ 
			sample_count = sample_size 
		}
		
		pop_proportion = BAF_level
		
		if( BAF_level < 0 ){ 
			pop_proportion = 0 
		}
		
		if( BAF_level > 1 ){ 
			pop_proportion = 1 
		}
		
		pval = calc_Pvalue_binomial_twotailed( sample_count, sample_size, pop_proportion )
		
		if( abs( BAF.mean - BAF_level ) < maxdist_BAF ) {
			pval=1
		}
		
	}else
	{
		pval = 0
		
	}	
	
	return( pval )
	
}

####################################################################################################
#' Calculate a p-value for a LogR value
#' TODO: this function is not used in Battenberg
#' @noRd
calc_LogR_Pvalue <-function( LogR, maxdist_LogR, LogR_level ) # kjd 27-2-2014
{
	if( is.finite( LogR_level ) )
	{
		pval = 0
		
		if( abs( LogR - LogR_level ) < maxdist_LogR ) {
			pval=1
		}
		
	}else
	{
		pval = 0
		
	}	
	
	return( pval )
	
}

#' Helper function to estimate rho from a given copy number state and it's BAF. The LogR is not used.
#' @noRd
estimate_rho <-function( LogR_value, BAFreq_value, nA_value, nB_value ) # kjd 10-3-2014
{
  rho_value = (2*BAFreq_value-1)/(2*BAFreq_value-BAFreq_value*(nA_value+nB_value)-1+nA_value)
  return( rho_value )
  
}

####################################################################################################
#' Helper function to calculate psi from a copy number fit, BAF, LogR, rho and a platform gamma
#' @noRd
estimate_psi <-function( LogR_value, BAFreq_value, nA_value, nB_value, rho_value, gamma_param ) # kjd 10-3-2014
{
  temp_value = 2^( - LogR_value / gamma_param )
  temp_value = temp_value * ( 2 + ( rho_value * ( nA_value + nB_value - 2 ) ) )
  #return(temp_value) # DCW this returns psi rather than psi_t, i.e. the average ploidy of normal and tumour cells
  temp_value = temp_value - ( 2 * ( 1 - rho_value ) )   
  psi_value = temp_value / rho_value    
  return( psi_value )	
}

#' Function that calculates rho and psi from a given reference segment, defined by ref_seg, with copy number state nA_ref and nB_ref
#' @noRd
get.psi.rho.from.ref.seg <-function( ref_seg, s, nA_ref, nB_ref, gamma_param = 1)
{
	BAFreq = s[ ref_seg, "b" ]
	LogR = s[ ref_seg, "r" ]
	
	rho = estimate_rho( LogR, BAFreq, nA_ref, nB_ref )
	psi = estimate_psi( LogR, BAFreq, nA_ref, nB_ref, rho, gamma_param )
		
	# ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
	nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
	nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
	ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"])
		
	# TODO DEBUG
	if (rho > 0) {
		ref_segment_info = list( psi = psi, rho = rho, ploidy = ploidy )
	} else {
		ref_segment_info = list( psi = NA, rho = NA, ploidy = NA )
	}

	
	
	return( ref_segment_info )	
}

#' This function decides if a segment is "clonal" (= TRUE) or not (= FALSE).
#' (The alternative hypothesis is that the tumour genome segment in question exhibits "sub-clonal" variation.)
#' We test the integer solutions for all 4 corners. Also, along side the hypothesis test for the BAF.
#' We use a decision rule based on LogR (we could use a hypothesis test which takes account of the variance in LogR, or a fixed “tolerance”).
#' If the null hypothesis is accepted for at least one corner, then we accept that
#' the tumour genome segment in question is "clonal".
#' @noRd
is.segment.clonal <-function( LogR, BAFreq, BAF.length, BAF.size, BAF.mean, BAF.sd, read_depth, rho, psi, gamma_param, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR ) # kjd 21-2-2014
{	
  # TODO: read_depth, siglevel_LogR and maxdist_LogR are no longer in use
  
	#270314 no longer used
	#pooled_BAF.size = read_depth * BAF.size
	
	# if we don't have a value for LogR, fill in 0
	if (is.na(LogR)) {
		LogR = 0
	}
  
  nA = (rho-1-(BAFreq-1)*2^(LogR/gamma_param)*((1-rho)*2+rho*psi))/rho
  nB = (rho-1+BAFreq*2^(LogR/gamma_param)*((1-rho)*2+rho*psi))/rho

  # if (any(is.na(nA) | is.na(nB)) | any(nA < 0 | nB < 0)) {
  #   # Reset any negative copy number to 0
  #   index = which(is.na(nA) | is.na(nB) | nA < 0 | nB < 0)
  #   print(paste("is.segment.clonal: Found negative copy number for segment", index, "BAF:", BAFreq[index], "logR:", LogR[index], "seg size:", BAF.size[index], "baf.sd:", BAF.sd[index]))
  #   nA[nA < 0 | is.na(nA)] = 0
  #   nB[nB < 0 | is.na(nB)] = 0
  # }


  nMajor = max(nA,nB, na.rm=T)
  nMinor = min(nA,nB, na.rm=T)

	# check for big shifts in nMajor - if there's a big shift, we shouldn't trust a clonal call
	nMajor.saved = nMajor
	## to make sure we're always in a positive square:
	#if(nMajor < 0) {
	#	nMajor = 0.01
	#}
	#
	#if(nMinor < 0) {
	#	nMinor = 0.01
	#}
	#DCW - increase nMajor and nMinor together, to avoid impossible combinations (with negative subclonal fractions)
	if(nMinor<0){
		if(BAFreq==1){
			#avoid calling infinite copy number
			nMajor = 1000
		}else{
			nMajor = nMajor + BAFreq * (0.01 - nMinor) / (1-BAFreq)
			if (nMajor<0) nMajor=1000
		}
		nMinor = 0.01		
	}
	
	# note that these are sorted in the order of ascending BAF:
	nMaj = c(floor(nMajor),ceiling(nMajor),floor(nMajor),ceiling(nMajor))
	nMin = c(ceiling(nMinor),ceiling(nMinor),floor(nMinor),floor(nMinor))
	x = floor(nMinor)
	y = floor(nMajor)
	
	# total copy number, to determine priority options
	ntot = nMajor + nMinor
	
	BAF_levels = (1-rho+rho*nMaj)/(2-2*rho+rho*(nMaj+nMin))
  #problem if rho=1 and nMaj=0 and nMin=0
  BAF_levels[nMaj==0 & nMin==0] = 0.5
	
	LogR_levels = gamma_param * log( (2-2*rho+rho*(nMaj+nMin))/(2-2*rho+rho*psi) , 2 ) # kjd 21-2-2014


  #DCW - just test corners on the nearest edge to determine clonality
  #If the segment is called as subclonal, this is the edge that will be used to determine the subclonal proportions that are reported first
  all.edges = orderEdges(BAF_levels, BAFreq, ntot,x,y)

  nMaj.test = all.edges[1,c(1,3)]
  nMin.test = all.edges[1,c(2,4)]
  test.BAF_levels = (1-rho+rho*nMaj.test)/(2-2*rho+rho*(nMaj.test+nMin.test))
  #problem if rho=1 and nMaj=0 and nMin=0
  test.BAF_levels[nMaj.test==0 & nMin.test==0] = 0.5
    
  whichclosestlevel.test = which.min(abs(test.BAF_levels-BAFreq))
  
  #270713 - problem caused by segments with constant BAF (usually 1 or 2)
  if(BAF.sd==0){
	  pval=0
  }else{
  	#pval[i] = t.test(BAFreq,alternative="two.sided",mu=BAF_levels[whichclosestlevel])$p.value
  	#pval = t.test(BAFreq,alternative="two.sided",mu=test.BAF_levels[whichclosestlevel.test])$p.value
  	pval = calc_Pvalue_t_twotailed( BAF.size, BAFreq, BAF.sd, test.BAF_levels[whichclosestlevel.test], maxdist_BAF)
  }
  #not necessary, because checked in calc_Pvalue_t_twotailed
  #if(min(abs(l-test.BAF_levels[whichclosestlevel.test]))<maxdist_BAF) {
  #  pval=1
  #}
  balanced = nMaj.test[whichclosestlevel.test] == nMin.test[whichclosestlevel.test]
  
  is.clonal = (pval > siglevel_BAF)
	# check for big shifts in nMajor - if there's a big shift, we shouldn't trust a clonal call
	# This is particularly problematic for very high cellularity samples, like some of the ovarian samples
	is.clonal = (pval > siglevel_BAF & nMajor - nMajor.saved <1)
  
  segment_info = list( is.clonal = is.clonal, balanced = balanced, nMaj.test = nMaj.test[whichclosestlevel.test] , nMin.test = nMin.test[whichclosestlevel.test] )
	
	return( segment_info )
	
}

####################################################################################################
#' This function calculates a t variate.
#' @noRd
calc_standardised_error <-function( LogR, BAFreq, BAF.length, BAF.size, BAF.mean, BAF.sd, rho, psi, gamma_param, maxdist_BAF ) # kjd 31-1-2014
{
	
	# if we don't have a value for LogR, fill in 0
	if (is.na(LogR)) {
		LogR = 0
	}
	nMajor = (rho-1+BAFreq*psi*2^(LogR/gamma_param))/rho
	nMinor = (rho-1+(1-BAFreq)*psi*2^(LogR/gamma_param))/rho
	
	# to make sure we're always in a positive square:
	if(nMajor < 0 | is.na(nMajor)) {
		nMajor = 0.01
	}
	
	if(nMinor < 0 | is.na(nMinor)) {
		nMinor = 0.01
	}
	
	# note that these are sorted in the order of ascending BAF:
	nMaj = c(floor(nMajor),ceiling(nMajor),floor(nMajor),ceiling(nMajor))
	nMin = c(ceiling(nMinor),ceiling(nMinor),floor(nMinor),floor(nMinor))
	x = floor(nMinor)
	y = floor(nMajor)
	
	# total copy number, to determine priority options
	ntot = nMajor + nMinor
	
	index_vect = which( (2-2*rho+rho*(nMaj+nMin)) != 0 )  # kjd 13-1-2014
	nMaj = nMaj[ index_vect ]  # kjd 13-1-2014
	nMin = nMin[ index_vect ]  # kjd 13-1-2014
	BAF_levels = (1-rho+rho*nMaj)/(2-2*rho+rho*(nMaj+nMin))
	
	whichclosestlevel = which.min(abs(BAF_levels-BAFreq))
	# if 0.5 and there are multiple options, finetune, because a random option got chosen 
	if( length( BAF_levels ) >= 3 ) {  # kjd 13-1-2014
		if (BAF_levels[whichclosestlevel]==0.5 && BAF_levels[2]==0.5 && BAF_levels[3]==0.5) {
			whichclosestlevel = ifelse(ntot>x+y+1,2,3)
		}  
	}  # kjd 13-1-2014  
	
	mu=BAF_levels[whichclosestlevel] # kjd 28-1-2014
	included_segment = 0 # kjd 31-1-2014
	if( BAF.size>0 ) { # kjd 13-1-2014

		if( BAF.sd==0 | length(mu)==0) {
			# pval=0 # kjd 31-1-2014
			tvar=0 # kjd 31-1-2014
			
		}else{
			# pval = t.test(BAFke,alternative="two.sided",mu=BAF_levels[whichclosestlevel])$p.value
			pval = calc_Pvalue_t_twotailed( BAF.size, BAF.mean, BAF.sd, mu, maxdist_BAF ) # kjd 31-1-2014
			
			tvar = studentise( BAF.size, BAF.mean, BAF.sd, mu ) # kjd 31-1-2014
			
			included_segment = 1 # kjd 31-1-2014
			
		}
	}else{ # kjd 13-1-2014
		# pval = 1 # kjd 13-1-2014 # kjd 31-1-2014
		tvar=0 # kjd 31-1-2014
		
	} # kjd 13-1-2014
		
	standard_error_info = list( included_segment = included_segment , tvar = tvar ) # kjd 31-1-2014
	
	return( standard_error_info )
	
}

####################################################################################################
#' This function computes various "distances", which are used as penalties for a copy number solution.
#' This function is called when searching for a clonal copy number solution.
#' One such distance is an estimate of the proportion of the tumour genome which is clonal.
#' For each segment of the genome, we test the null hypothesis is that
#' the tumour genome segment in question is "clonal". The alternative hypothesis is that
#' the tumour genome segment in question exhibits "sub-clonal" variation.	  
#' @noRd
calc_distance <-function( segs, dist_choice, rho, psi, gamma_param, uninformative_BAF_threshold=0.51 ) # kjd 10-2-2014  
{
  s = segs
  
  if( dist_choice == 0 ) # original ASCAT distance
  {
    nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
    nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
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
    #dist_value = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"], na.rm=T)
    #DCW 310314 - retry weighting
    dist_value = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"] * ifelse(s[,"b"]<=uninformative_BAF_threshold,0.05,1), na.rm=T)
    
    minimise = TRUE
    
  }else if( dist_choice == 1 ){ # new similarity measure suggested by DW 7-3-2014
    nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
    nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
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
    # dist_value = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"], na.rm=T)
    
    dist_value = sum((0.5-abs(nMinor - pmax(round(nMinor),0)))^2 * s[,"length"], na.rm=T)
    
    minimise = FALSE
    
  }else if( dist_choice == 2 ){ # adapted DW's 7-3-2014 measure by SD 8-8-2014 that takes into account both major and minor alleles
    nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
    nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
    # choose the minor allele
    nMinor = NULL
    nMajor = NULL
    if (sum(nA,na.rm=T) < sum(nB,na.rm=T)) {
      nMinor = nA
      nMajor = nB
    }
    else {
      nMinor = nB
      nMajor = nA
    }
    
    dist_value = 0.5*sum(((0.5-abs(nMinor - pmax(round(nMinor),0)))^2 + (0.5-abs(nMajor - pmax(round(nMajor),0)))^2) * s[,"length"], na.rm=T)
    
    minimise = FALSE
    
  }else if( dist_choice == 3 ){ # adapted DW's 7-3-2014 measure by SD 8-8-2014 that takes into account both major and minor alleles and takes the mean, while it also penalises for the number of homozygous deletions
    nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
    nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
    # choose the minor allele
    nMinor = NULL
    nMajor = NULL
    if (sum(nA,na.rm=T) < sum(nB,na.rm=T)) {
      nMinor = nA
      nMajor = nB
    }
    else {
      nMinor = nB
      nMajor = nA
    }
    
    # Penalise homozygous deletions twice as hard as other segments
    # - the penalty term is increased to make it less likely that hom dels occur
    # - the segment length is increased to penalise harder for longer segments
    segs_penalty = (0.5-abs(nMinor - pmax(round(nMinor),0)))^2 + (0.5-abs(nMajor - pmax(round(nMajor),0)))^2
    hom_del = nMinor<0.5 & nMajor<0.5 & nMinor>=0 & nMajor>=0
    segs_penalty[which(hom_del)] = segs_penalty[which(hom_del)]*4
    
    dist_value = 0.5*sum(segs_penalty * (s[,"length"] * ifelse(hom_del, 2, 1)), na.rm=T)
    
    minimise = FALSE
  }
  
  distance_info = list( distance_value = dist_value , minimise = minimise )
  
  return( distance_info )
}

####################################################################################################
#' This function computes various "distances", which are used as penalties for a copy number solution
#' One such distance is an estimate of the proportion of the tumour genome which is clonal.
#' For each segment of the genome, we test the null hypothesis is that
#' the tumour genome segment in question is "clonal". The alternative hypothesis is that
#' the tumour genome segment in question exhibits "sub-clonal" variation.    
#' @noRd
calc_distance_clonal <-function( segs, dist_choice, rho, psi, gamma_param, read_depth, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR, uninformative_BAF_threshold) # kjd 10-2-2014  
{
  s = segs
  
  pval = NULL
  
  # BAFpvals = vector(length=length(BAFseg))
  
  genome_size = 0
  clonal_genome_size = 0
  seg_count = 0 # kjd 24-1-2014
  clonal_seg_count = 0 # kjd 24-1-2014
  n_included_segments = 0 # kjd 31-1-2014
  included_genome_size = 0 # kjd 31-1-2014
  sum1 = 0 # kjd 31-1-2014
  sum2 = 0 # kjd 31-1-2014
  sum3 = 0 # kjd 31-1-2014
  sum_ln_lratio = 0 # kjd 10-2-2014
  
  max_clonal_segment = 0 # There may be no clonal segments, in which case this remains zero.
  max_clonal_segment_size = 0
  
  ref_maj = NA
  ref_min = NA
  
  for(i in 1:nrow(s)) {
    
    BAFreq = s[ i, "b" ] # l = BAFlevels[i]
    
    if( BAFreq > uninformative_BAF_threshold )
    {
      LogR = s[ i, "r" ]
      
      BAF.length = s[ i, "length" ]
      BAF.size = s[ i, "size" ]
      BAF.mean = s[ i, "mean" ]
      BAF.sd = s[ i, "sd" ]
      
      #
      # Calculate P values
      #
      
      segment_info = is.segment.clonal( LogR, BAFreq, BAF.length, BAF.size, BAF.mean, BAF.sd, read_depth, rho, psi, gamma_param, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR ) # kjd 21-2-2014
      is.clonal = segment_info$is.clonal # kjd 21-2-2014			
      
      nMaj = segment_info$nMaj
      nMin = segment_info$nMin			
      is.balanced = segment_info$balanced
      
      segment_size = BAF.length # OR segment_size = BAF.size ?
      genome_size = genome_size + segment_size
      seg_count = seg_count + 1 # kjd 24-1-2014
      
      # if( pval[i] > siglevel_BAF ){
      if(is.clonal){ # kjd 21-2-2014
        clonal_genome_size = clonal_genome_size + segment_size
        clonal_seg_count = clonal_seg_count + 1 # kjd 24-1-2014
        
        if( max_clonal_segment_size < segment_size & !is.balanced) #balanced check added by DCW 160314
        {
          max_clonal_segment = i
          max_clonal_segment_size = segment_size
          
          ref_maj = nMaj
          ref_min = nMin
        }
        
      }
      
      #
      # Calculate "standardised error"
      #
      
      standard_error_info = calc_standardised_error( LogR, BAFreq, BAF.length, BAF.size, BAF.mean, BAF.sd, rho, psi, gamma_param, maxdist_BAF ) # kjd 31-1-2014
      
      included_segment = standard_error_info$included_segment # kjd 31-1-2014
      tvar = standard_error_info$tvar # kjd 31-1-2014
      
      n_included_segments = n_included_segments + included_segment # kjd 31-1-2014
      if( included_segment > 0 )
      {
        included_genome_size = included_genome_size + segment_size # kjd 31-1-2014
        
      }
      sum1 = sum1 + tvar^2 # kjd 31-1-2014
      
      sum2 = sum2 + ( BAFreq - BAF.mean )^2
      
      sum3 = sum3 + ( segment_size * ( BAFreq - BAF.mean )^2 )
      
      #
      # Calculate log likelihood ratio
      #
      
      ln_lratio = calc_ln_likelihood_ratio( LogR, BAFreq, BAF.length, BAF.size, BAF.mean, read_depth, rho, psi, gamma_param, maxdist_BAF ) # kjd 10-2-2014
      
      sum_ln_lratio = sum_ln_lratio + ln_lratio
      
    }
    
  }
  
  #
  # Calculate proportion of genome which is "clonal":
  #
  
  clonal_proportion = 0
  if( genome_size > 0 ){
    clonal_proportion = clonal_genome_size / genome_size
    
  }
  
  #
  # Calculate "distances":
  #
  
  dist1 = 0 # kjd 3-2-2014
  if( n_included_segments > 0 ){
    dist1 = sum1 / n_included_segments
    
  } # kjd 3-2-2014
  
  dist2 = 0 # kjd 3-2-2014
  if( seg_count > 0 ){ 
    dist2 = sum2 / seg_count
    
  } # kjd 3-2-2014
  
  dist3 = 0  # kjd 3-2-2014
  if( genome_size > 0 ){
    dist3 = sum3 / genome_size
    
  } # kjd 3-2-2014
  
  
  
  if( dist_choice == 0 )
  {
    dist_value = clonal_proportion
    minimise = FALSE
  }
  
  if( dist_choice == 1 )
  {
    dist_value = dist1
    minimise = TRUE
  }
  
  if( dist_choice == 2 )
  {
    dist_value = dist2
    minimise = TRUE
  }
  
  if( dist_choice == 3 )
  {
    dist_value = dist3
    minimise = TRUE
  }
  
  if( dist_choice == 4 )
  {
    dist_value = sum_ln_lratio
    minimise = FALSE
  }
  
  distance_info = list( distance_value = dist_value , minimise = minimise , max_clonal_segment = max_clonal_segment, ref_maj = ref_maj, ref_min = ref_min ) # kjd 10-2-2014
  
  # return( clonal_proportion ) # kjd 24-1-2014
  
  return( distance_info ) # kjd 10-2-2014
  
}

#' Function extends the ASCAT \code{make_segments} function to make segments
#' of constant BAF and LogR. This function returns a matrix with for each
#' segment the LogR, BAF, the length of the segment (twice), and the mean and 
#' standard deviation of the BAF values
#' @noRd
get_segment_info = function(segLogR , segBAF.table) {
  segBAF = segBAF.table[,5]
  
  names(segBAF) = rownames(segBAF.table)
  names(segLogR) = rownames(segBAF.table)
  
  b = segBAF
  r = segLogR[names(segBAF)]
  pcf_segments = ASCAT::make_segments(r,b)
  
#   m = matrix(ncol = 2, nrow = length(b))
#   m[,1] = r
#   m[,2] = b
#   m = as.matrix(na.omit(m))
#   pcf_segments = matrix(ncol = 3, nrow = dim(m)[1])
#   colnames(pcf_segments) = c("r","b","length");
#   index = 0;
#   previousb = -1;
#   previousr = 1E10;
#   for (i in 1:dim(m)[1]) {
#     if (m[i,2] != previousb || m[i,1] != previousr) {
#       index=index+1;
#       count=1;
#       pcf_segments[index, "r"] = m[i,1];
#       pcf_segments[index, "b"] = m[i,2];
#     }
#     else {
#       count = count + 1;
#     }
#     pcf_segments[index, "length"] = count;
#     previousb = m[i,2];
#     previousr = m[i,1];
#   }
#   
#   # pcf_segments = as.matrix(na.omit(pcf_segments))[,] # kjd 10-1-2014 This version caused bug in R on laptop.
#   pcf_segments = as.matrix(na.omit(pcf_segments)) # kjd 10-1-2014 This version resolved bug in R on laptop. (Problem with installed version of R?)
#   
  segs = matrix(ncol = 6, nrow = nrow(pcf_segments))
  colnames(segs) = c("r","b","length","size", "mean", "sd")
  segs[ , c("r","b","length")] = pcf_segments

  for( i in 1:nrow(segs) ) {
		BAFreq = segs[i, "b"] # l = BAFlevels[i]		
		index_vect = which( segBAF.table[ , 5] == BAFreq )
		BAFke = segBAF.table[index_vect, 4] # column 4 contains "phased BAF" values; # kjd 6-1-2014 
		
		segs[i, "size"] = length(BAFke)
		segs[i, "mean"] = mean(BAFke)
		segs[i, "sd"] = sd(BAFke)
  }
  return(segs);
}

####################################################################################################
#' Helper function to find new rho and psi boundaries given a current optimum pair. 
#' @noRd
get_new_bounds = function( input_optimum_pair, ininitial_bounds ) # kjd 21-2-2014
{
  
  psi_optimum = input_optimum_pair$psi
  rho_optimum = input_optimum_pair$rho
  
  psi_min_initial = ininitial_bounds$psi_min
  psi_max_initial = ininitial_bounds$psi_max
  rho_min_initial = ininitial_bounds$rho_min
  rho_max_initial = ininitial_bounds$rho_max
  
  psi_range = 0.1 * ( psi_max_initial - psi_min_initial )
  #rho_range = 0.1 * ( rho_max_initial - rho_min_initial )
  #DCW 170314 - rho range depends on optimum value of rho
  rho_range = 0.1 * rho_optimum
  
  if( (psi_optimum - 0.5 * psi_range) < psi_min_initial )
  {
    psi_min = psi_min_initial
    psi_max = psi_min_initial + psi_range
    
  }else
  {
    if( (psi_optimum + 0.5 * psi_range) > psi_max_initial )
    {
      psi_min = psi_max_initial - psi_range
      psi_max = psi_max_initial
      
    }else
    {
      psi_min = psi_optimum - 0.5 * psi_range
      psi_max = psi_optimum + 0.5 * psi_range
      
    }
  }
  
  if( (rho_optimum - 0.5 * rho_range) < rho_min_initial )
  {
    rho_min = rho_min_initial
    rho_max = rho_min_initial + rho_range
    
  }else
  {
    if( (rho_optimum + 0.5 * rho_range) > rho_max_initial )
    {
      rho_min = rho_max_initial - rho_range
      rho_max = rho_max_initial
      
    }else
    {
      rho_min = rho_optimum - 0.5 * rho_range
      rho_max = rho_optimum + 0.5 * rho_range
      
    }
  }
  
  new_bounds = list( psi_min = psi_min, psi_max = psi_max, rho_min = rho_min, rho_max = rho_max )
  
  
  return( new_bounds )
  
}
 
####################################################################################################
#' function to create the distance matrix (distance for a range of ploidy and tumor percentage values)
#' input: segmented LRR and BAF and the value for gamma_param
#' @noRd
create_distance_matrix = function(s, dist_choice, gamma_param, uninformative_BAF_threshold=0.51, min_rho=0.1, max_rho=1, min_psi=1, max_psi=5.4) {
  psi_pos = seq(min_psi,max_psi,0.05) 
  rho_pos = seq(min_rho,max_rho,0.01)
  d = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) = psi_pos
  colnames(d) = rho_pos
  dmin = 1E20;
  for(i in 1:length(psi_pos)) {
    psi = psi_pos[i]
    for(j in 1:length(rho_pos)) {
      rho = rho_pos[j]
      
      distance_info = calc_distance( s, dist_choice, rho, psi, gamma_param, uninformative_BAF_threshold=uninformative_BAF_threshold ) # kjd 10-2-2014
	  
	  d[i,j] = distance_info$distance_value
      # minimise = distance_info$minimise
	  
    }
  }
  
  minimise = distance_info$minimise
  
  distance_matrix_info = list( distance_matrix = d , minimise = minimise )
  
  # return(d)
  return( distance_matrix_info )
  
}

#' Helper function to create the clonal distance matrix for a range of
#' rho and psi values
#' @noRd
create_distance_matrix_clonal = function( segs, dist_choice, gamma_param, read_depth, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR, uninformative_BAF_threshold, new_bounds) # kjd 18-12-2013
{
  psi_min = new_bounds$psi_min
  psi_max = new_bounds$psi_max
  rho_min = new_bounds$rho_min
  rho_max = new_bounds$rho_max
  
  s = segs
    
  psi_range = psi_max - psi_min
  rho_range = rho_max - rho_min
  
  delta_psi = psi_range / 100
  delta_rho = rho_range / 100
  
  psi_pos = seq( psi_min, psi_max, delta_psi )
  rho_pos = seq( rho_min, rho_max, delta_rho )
  
  # psi_pos = seq(1,5.4,0.05) 
  # rho_pos = seq(0.1,1.05,0.01)
  
  ref_seg_matrix = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  ref_major = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  ref_minor = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(ref_seg_matrix) = psi_pos
  colnames(ref_seg_matrix) = rho_pos
  rownames(ref_major) = psi_pos
  colnames(ref_major) = rho_pos
  rownames(ref_minor) = psi_pos
  colnames(ref_minor) = rho_pos
  
  d = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) = psi_pos
  colnames(d) = rho_pos
  # dmin = 1E20;
  for(i in 1:length(psi_pos)) {
    psi = psi_pos[i]
    for(j in 1:length(rho_pos)) {
      rho = rho_pos[j]      
            
      # clonal_proportion = calc_clonal_proportion( s, LogRvals, BAFvals, segBAF.table, rho, psi, gamma_param, siglevel_BAF, maxdist_BAF ) # kjd 18-12-2013
	  distance_info = calc_distance_clonal( s, dist_choice, rho, psi, gamma_param, read_depth, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR, uninformative_BAF_threshold) # kjd 10-2-2014
	  
	  distance_value = distance_info$distance_value # kjd 10-2-2014
      # minimise = distance_info$minimise # kjd 10-2-2014
	  max_clonal_segment = distance_info$max_clonal_segment
	  
      d[i,j] = distance_value # kjd 10-2-2014
      ref_seg_matrix[i,j] = max_clonal_segment
 	
 	  ref_major[i,j] = distance_info$ref_maj
 	  ref_minor[i,j] = distance_info$ref_min
    }
  }
  
  minimise = distance_info$minimise # kjd 10-2-2014
  
  distance_matrix_info = list( distance_matrix = d , minimise = minimise , ref_seg_matrix = ref_seg_matrix, ref_major = ref_major, ref_minor = ref_minor ) # kjd 10-2-2014
  
  # return(d) # kjd 10-2-2014
  return( distance_matrix_info ) # kjd 10-2-2014
  
}

####################################################################################################
#' Helper function to calculate a square distance
#' @noRd
calc_square_distance <-function( pt1, pt2 ) # kjd 27-2-2014
{
	dsqr = ( pt1[1] - pt2[1] )^2 + ( pt1[2] - pt2[2] )^2
	
	return( dsqr )
	
}

####################################################################################################
#' This function is an alternative procedure for finding the optimum (psi, rho) pair.
#' This function first finds all the find all the global optima,
#' and then finds the centroid of this set of globla optima.
#' Then we find the global optimum which is nearest to the centroid. 
#' (When the set of global optima is convex, we expect the selected optimum to be at the centroid.)
#' @param d A distance matrix
#' @param ref_seg_matrix The corresponding ref seg matrix that belongs to d
#' @param ref_major The corresponding major allele values with d
#' @param ref_minor The corresponding minor allele values with d
#' @param s A segmented BAF/LogR data.frame from \code{get_segment_info}
#' @param dist_choice Some distance metrics require adaptation of the data (i.e. log transform) 
#' @param minimise Boolean whether we're minimising or maximising
#' @param new_bounds The rho/psi boundaries between we are searching for a solution. This is a named list with values psi_min, psi_max, rho_min, rho_max
#' @param distancepng String where the sunrise distance plot will be saved
#' @param gamma_param The platform gamma
#' @param siglevel_BAF The level at which BAF becomes significant TODO: this option is no longer used
#' @param maxdist_BAF TODO: this option is no longer used
#' @param siglevel_LogR The p-value at which logR becomes significant when establishing whether a segment should be subclonal
#' @param maxdist_LogR The maximum distance allowed as slack when establishing the significance. This allows for the case when a breakpoint is missed, the segment would then not automatically become subclonal
#' @param allow100percent Boolean whether to allow for a 100"\%" cellularity solution
#' @param uninformative_BAF_threshold The threshold above which BAF becomes uninformative
#' @param read_depth TODO: this option is no longer used
#' @return A list with fields optima_info_without_ref and optima_info
#' @export
find_centroid_of_global_minima <- function( d, ref_seg_matrix, ref_major, ref_minor, s, dist_choice, minimise, new_bounds, distancepng, gamma_param, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR, allow100percent, uninformative_BAF_threshold, read_depth) # kjd 28-2-2014
{

  #Theoretmaxdist_BAF = sum(rep(0.25,dim(s)[1]) * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1),na.rm=T)
  #DCW 180711 - try weighting BAF=0.5 equally with other points
  # Theoretmaxdist_BAF = sum(rep(0.25,dim(s)[1]) * s[,"length"],na.rm=T)
  
  
  if( !(minimise) ) # kjd 12-2-2013
  {
	d = - d # This ensures that we "maximise" instead of "minimise"!	
  }  
  
  # Find height of global minima;
  # (subject to additional conditions: percentzero > 0.01 | perczeroAbb > 0.1)
  
  gmin = max( d )
  for (i in 1:(dim(d)[1])) {
    for (j in 1:(dim(d)[2])) {
      psi = as.numeric(rownames(d)[i])
      rho = as.numeric(colnames(d)[j])
      nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
      nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
      
      # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
      ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
      
      percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
      perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
      # the next can happen if BAF is a flat line at 0.5 
      if (is.na(perczeroAbb)) {
          perczeroAbb = 0
      }
      
      # commented out by kjd 6-3-2014
       #if( percentzero > 0.01 | perczeroAbb > 0.1 ) { # kjd 6-3-2014
          
          if( d[i,j] <= gmin ) {
              gmin = d[i,j]
              
          }         
       #}            
    }
  }
  
  # Find all global minima;
  # (subject to additional conditions: percentzero > 0.01 | perczeroAbb > 0.1)
  
  nropt = 0
  localmin = NULL
  optima = list()
  
  for (i in 1:(dim(d)[1])) {
    for (j in 1:(dim(d)[2])) {
      if( d[i,j] == gmin ) {
        psi = as.numeric(rownames(d)[i])
        rho = as.numeric(colnames(d)[j])
        nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
        nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma_param)*((1-rho)*2+rho*psi))/rho
        
        # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
        ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);
      
        percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
        perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
        # the next can happen if BAF is a flat line at 0.5 
        if (is.na(perczeroAbb)) {
          perczeroAbb = 0
        }

        # goodnessOfFit = (1-m/Theoretmaxdist_BAF) * 100
        goodnessOfFit = gmin #DCW 250314 goodnessOfFit is the same as gmin, because the metric is the total amount of the genome that is clonal
        nropt = nropt + 1
        optima[[nropt]] = c(gmin,i,j,ploidy,goodnessOfFit)
        localmin[nropt] = gmin
        
      }     
    }
  }	  

  #
  # Find a "centroid" of the set of global minima:
  #
  
  grid_x_vect = unlist( lapply( optima , function(z){ z[2] } ) )
  grid_y_vect = unlist( lapply( optima , function(z){ z[3] } ) )
  
  centre_x = mean( median( grid_x_vect ) )
  centre_y = mean( median( grid_y_vect ) )
  
  centre = c( centre_x, centre_y )
      	
	index = 1
	sqrdist_min = (dim(d)[1])^2 + (dim(d)[2])^2
	for (i in 1:length(optima)) {
		
		grid_x = optima[[i]][2] # grid_i = ( psi_opt1 - 1 ) * 20
		grid_y = optima[[i]][3] # grid_j = ( rho_opt1 - 0.1 ) * 100
		
		grid_point = c( grid_x, grid_y )
		
		sqrdist = calc_square_distance( grid_point, centre )
	    
		if( sqrdist <= sqrdist_min ) {
			sqrdist_min = sqrdist
			index = i
		}
    }	
	
	grid_x = optima[[index]][2] # grid_i = ( psi_opt1 - 1 ) * 20
	grid_y = optima[[index]][3] # grid_j = ( rho_opt1 - 0.1 ) * 100
	    
    psi_opt1 = as.numeric(rownames(d)[optima[[index]][2]])
    rho_opt1 = as.numeric(colnames(d)[optima[[index]][3]])
    if(rho_opt1 > 1) {
       rho_opt1 = 1
    }
    ploidy_opt1 = optima[[index]][4]
    goodnessOfFit_opt1 = optima[[index]][5]
    
    ref_seg = ref_seg_matrix[ grid_x, grid_y ]
    
    # store optima for plotting later
    rhos = rho_opt1
    psis = psi_opt1
      #
	  # Write to clonal info file:
	  #
	  
     if( isTRUE(minimise) ) # kjd 12-2-2013
     {
		dist_optima = gmin # when we "minimise";
		
	 }else
	 {
		dist_optima = - gmin # Recall that when we "maximise", we replace "d" by "-d";
		goodnessOfFit_opt1 = -goodnessOfFit_opt1 #DCW 250314		
	 }
	  
	  print(paste("goodnessOfFit from grid=",goodnessOfFit_opt1,sep=""))
	#DCW 140314
	optima_info_without_ref = list( nropt = nropt, psi_opt1 = psi_opt1, rho_opt1 = rho_opt1, ploidy_opt1 = ploidy_opt1, ref_seg = ref_seg, goodnessOfFit_opt1 = goodnessOfFit_opt1 )
	
	#DCW if no ref segment found, there is no tumour present
	if(ref_seg==0){
		psi_opt1 = 2
		rho_opt1 = 1
		ploidy_opt1=2
		goodnessOfFit_opt1 = 1
	}else{
		ref_segment_info = get.psi.rho.from.ref.seg( ref_seg, s, ref_major[ grid_x, grid_y ], ref_minor[ grid_x, grid_y ], gamma_param)
		
		psi_opt1 = ref_segment_info$psi
		rho_opt1 = ref_segment_info$rho
		ploidy_opt1 = ref_segment_info$ploidy
		
		# TODO DEBUG
		if (!is.na(rho_opt1)) {
		#goodness of fit is the same as the distance measure for fraction of genome that is clonal
  		distance.info = calc_distance_clonal( s, dist_choice, rho_opt1, psi_opt1, gamma_param, read_depth, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR, uninformative_BAF_threshold)
  		goodnessOfFit_opt1 = distance.info$distance_value		
		#goodnessOfFit_opt1 = ref_segment_info$goodnessOfFit_opt1
		} else {
			goodnessOfFit_opt1 = Inf
		}

	} 

	# store optima for plotting later
    	rhos = c(rhos, rho_opt1)
    	psis = c(psis, psi_opt1)
  
  # separated plotting from logic: create distanceplot here
  if (!is.na(distancepng)) {
    png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
  }
  clonal_findcentroid.plot(minimise, dist_choice, -d, psis, rhos, new_bounds)
  if (!is.na(distancepng)) { dev.off() }
	
	optima_info = list( nropt = nropt, psi_opt1 = psi_opt1, rho_opt1 = rho_opt1, ploidy_opt1 = ploidy_opt1, ref_seg = ref_seg, goodnessOfFit_opt1 = goodnessOfFit_opt1 ) # kjd 10-3-2014
	
	return( list(optima_info_without_ref=optima_info_without_ref, optima_info=optima_info) )
}

#' A modified ASCAT main function to fit Battenberg
#' 
#' This function returns an initial rho and psi estimate for a clonal copy number fit. It uses an internal distance metric to create a distance matrix.
#' Using that matrix it will search for a rho and psi combination that yields the least heavy penalty.
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param baf (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
#' @param lrrsegmented log R, segmented, in genomic sequence (all probes), with probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param chromosomes a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#' @param dist_choice The distance metric to be used internally to penalise a copy number solution
#' @param distancepng if NA: distance is plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param copynumberprofilespng if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param nonroundedprofilepng if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file (Default NA)
#' @param cnaStatusFile File where the copy number profile status is written to. This contains either the message "No suitable copy number solution found" or "X copy number solutions found" (Default copynumber_solution_status.txt)
#' @param gamma technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 "\%" aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays) (Default 0.55)
#' @param allow100percent A boolean whether to allow a 100"\%" cellularity solution
#' @param reliabilityFile String to where fit reliabilty information should be written. This file contains backtransformed BAF and LogR values for segments using the fitted copy number profile (Default NA)
#' @param min.ploidy The minimum ploidy to consider (Default 1.6)
#' @param max.ploidy The maximum ploidy to consider (Default 4.8)
#' @param min.rho The minimum cellularity to consider (Default 0.1)
#' @param max.rho The maximum cellularity to consider (Default 1.0)
#' @param min.goodness The minimum goodness of fit for a solution to have to be considered (Default 63)
#' @param uninformative_BAF_threshold The threshold beyond which BAF becomes uninformative (Default 0.51)
#' @param chr.names A vector with chromosome names used for plotting
#' @param analysis A String representing the type of analysis to be run, this determines whether the distance figure is produced (Default paired)
#' @return A list with fields psi, rho and ploidy
#' @export
#the limit on rho is lenient and may lead to spurious solutions
runASCAT = function(lrr, baf, lrrsegmented, bafsegmented, chromosomes, dist_choice, distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA, cnaStatusFile = "copynumber_solution_status.txt", gamma = 0.55, allow100percent,reliabilityFile=NA,min.ploidy=1.6,max.ploidy=4.8,min.rho=0.1,max.rho=1.0,min.goodness=63, uninformative_BAF_threshold = 0.51, chr.names, analysis="paired") {
  ch = chromosomes
  b = bafsegmented
  r = lrrsegmented[names(bafsegmented)]

  # Adapt the rho/psi boundaries for the local maximum searching below to work
  dist_min_psi = max(min.ploidy-0.6, 0)
  dist_max_psi = max.ploidy+0.6 
  dist_min_rho = max(min.rho-0.03, 0.05)
  dist_max_rho = max.rho+0.03
  
  s = ASCAT::make_segments(r,b)
  dist_matrix_info <- create_distance_matrix( s, dist_choice, gamma, uninformative_BAF_threshold=uninformative_BAF_threshold, min_psi=dist_min_psi, max_psi=dist_max_psi, min_rho=dist_min_rho, max_rho=dist_max_rho)  
  d = dist_matrix_info$distance_matrix
  minimise = dist_matrix_info$minimise

  #TheoretMaxdist = sum(rep(0.25,dim(s)[1]) * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1),na.rm=T)
  #DCW 180711 - try weighting BAF=0.5 equally with other points
  TheoretMaxdist = sum(rep(0.25,dim(s)[1]) * s[,"length"],na.rm=T)

  if( !(minimise) ) # kjd 10-3-2014
  {
	d = - d # This ensures that we "maximise" instead of "minimise"!	
  }
  
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
	ploidy_opt1 = ploidy
      
        percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
        perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
        # the next can happen if BAF is a flat line at 0.5 
        if (is.na(perczeroAbb)) {
          perczeroAbb = 0
        }

        #goodnessOfFit = (1-m/TheoretMaxdist) * 100
        #140314 - DCW
        if(minimise){
			goodnessOfFit = (1-m/TheoretMaxdist) * 100
		}else{
			goodnessOfFit = -m/TheoretMaxdist * 100 # we have to use minus to reverse d=-d above
		}
		
		print(paste("ploidy=",ploidy,",rho=",rho,",goodness=",goodnessOfFit,",percentzero=",percentzero,", perczerAbb=",perczeroAbb,sep=""))
		if (ploidy >= min.ploidy & ploidy <= max.ploidy & rho >= min.rho & goodnessOfFit >= min.goodness & (percentzero > 0.01 | perczeroAbb > 0.1)) {	
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

          #goodnessOfFit = (1-m/TheoretMaxdist) * 100
        #140314 - DCW
        if(minimise){
			goodnessOfFit = (1-m/TheoretMaxdist) * 100
		}else{
			goodnessOfFit = -m/TheoretMaxdist * 100 # we have to use minus to reverse d=-d above
		}

          if (ploidy > min.ploidy & ploidy < max.ploidy & rho >= min.rho & goodnessOfFit >= min.goodness) {
            nropt = nropt + 1
            optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
            localmin[nropt] = m
          }
        }
      }
    }
  }

  # added for output to plotting
  psi_opt1_plot = vector(mode="numeric")
  rho_opt1_plot = vector(mode="numeric")

  if (nropt>0) {
    write.table(paste(nropt, " copy number solutions found", sep=""), file=cnaStatusFile, quote=F, col.names=F, row.names=F)
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
        psi_opt1_plot = c(psi_opt1_plot, psi_opt1)
        rho_opt1_plot = c(rho_opt1_plot, rho_opt1)
        # points((psi_opt1-1)/4.4,(rho_opt1-0.1)/0.95,col="green",pch="X", cex = 2)
      }
    }
  } else {
    write.table(paste("no copy number solutions found", sep=""), file=cnaStatusFile, quote=F, col.names=F, row.names=F)
    print("No suitable copy number solution found")
    psi = NA
    ploidy = NA
    rho = NA
    psi_opt1_plot = -1
    rho_opt1_plot = -1
  }
  
  # NAP: only create this plot for 'paired' analysis mode and not cell_line or germline; it shows strange behaviour and halts execution
  if (analysis=="paired"){	
  # separated plotting from logic: create distanceplot here
  if (!is.na(distancepng)) {
    png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
  }
  ASCAT::ascat.plotSunrise(-d, psi_opt1_plot, rho_opt1_plot,minimise)
  if (!is.na(distancepng)) { dev.off() }
}

  if(nropt>0) {

    rho = rho_opt1
    psi = psi_opt1
    ploidy = ploidy_opt1
    
    nAfull = (rho-1-(b-1)*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nBfull = (rho-1+b*2^(r/gamma)*((1-rho)*2+rho*psi))/rho
    nA = pmax(round(nAfull),0)
    nB = pmax(round(nBfull),0)
  
    rBacktransform = gamma*log((rho*(nA+nB)+(1-rho)*2)/((1-rho)*2+rho*psi),2)
    bBacktransform = (1-rho+rho*nB)/(2-2*rho+rho*(nA+nB))
    rConf = ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
    bConf = ifelse(bBacktransform!=0.5,pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))),NA)
    #DCW 150711 - get deviations from expected values
    if(!is.na(reliabilityFile)){
		write.table(data.frame(segmentedBAF=b,backTransformedBAF=bBacktransform,confidenceBAF=bConf,segmentedR=r,backTransformedR=rBacktransform,confidenceR=rConf,nA=nA,nB=nB,nAfull=nAfull,nBfull=nBfull), reliabilityFile,sep=",",row.names=F)
	}
    confidence = ifelse(is.na(rConf),bConf,ifelse(is.na(bConf),rConf,(rConf+bConf)/2))    

	# Create plot
	if (!is.na(copynumberprofilespng)) {
	  png(filename = copynumberprofilespng, width = 2000, height = 500, res = 200)
	}
	ASCAT::ascat.plotAscatProfile(n1all = nA, n2all = nB, heteroprobes = TRUE, ploidy = ploidy_opt1, rho = rho_opt1, goodnessOfFit = goodnessOfFit_opt1, nonaberrant = FALSE, ch = ch, lrr = lrr, bafsegmented = bafsegmented, chrs=chr.names)
	if (!is.na(copynumberprofilespng)) { dev.off() }
	
	# separated plotting from logic: create nonrounded copy number profile plot here
	if (!is.na(nonroundedprofilepng)) {
	  png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200)
	}
	# clonal_runascat.plot3(rho_opt1, goodnessOfFit_opt1, ploidy_opt1, nAfull, nBfull, ch, lrr, bafsegmented)
	ASCAT::ascat.plotNonRounded(ploidy = ploidy_opt1, rho = rho_opt1, goodnessOfFit = goodnessOfFit_opt1, nonaberrant = FALSE, nAfull = nAfull, nBfull = nBfull, bafsegmented = bafsegmented, ch = ch, lrr = lrr, chrs=chr.names)
	if (!is.na(nonroundedprofilepng)) { dev.off() }
  
  }
  output_optimum_pair = list(psi = psi, rho = rho, ploidy = ploidy)
  return( output_optimum_pair ) # kjd 20-2-2014 
}

####################################################################################################
#' ASCAT like function to obtain a clonal copy number profile
#' 
#' This function takes an initial optimum rho/psi pair and uses
#' an internal distance metric to calculate a score for each rho/psi pair allowed.
#' The solution with the best score is then taken to obtain a global copy number
#' profile. This function performs both a grid search and tries to find a reference
#' segment, but the grid search result is always used for now.
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param baf (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
#' @param lrrsegmented log R, segmented, in genomic sequence (all probes), with probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param chromosomes a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#' @param segBAF.table Segmented BAF data.frame from \code{get_segment_info}
#' @param input_optimum_pair A list containing fields for rho, psi and ploidy, as is output from \code{runASCAT}
#' @param dist_choice The distance metric to be used internally to penalise a copy number solution
#' @param distancepng if NA: distance is plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param copynumberprofilespng if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param nonroundedprofilepng if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file (Default NA)
#' @param gamma_param technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 "\%" aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays) (Default 0.55)
#' @param read_depth TODO: unused parameter that should be removed
#' @param uninformative_BAF_threshold The threshold beyond which BAF becomes uninformative
#' @param allow100percent A boolean whether to allow a 100"\%" cellularity solution
#' @param reliabilityFile String to where fit reliabilty information should be written. This file contains backtransformed BAF and LogR values for segments using the fitted copy number profile (Default NA)
#' @param psi_min_initial Minimum psi value to be considered (Default: 1.0)
#' @param psi_max_initial Maximum psi value to be considered (Default: 5.4)
#' @param rho_min_initial Minimum rho value to be considered (Default: 0.1)
#' @param rho_max_initial Maximum rho value to be considered (Default: 1.05)
#' @param chr.names A vector with chromosome names used for plotting
#' @return A list with fields output_optimum_pair, output_optimum_pair_without_ref, distance, distance_without_ref, minimise and is.ref.better
#' @export
run_clonal_ASCAT = function(lrr, baf, lrrsegmented, bafsegmented, chromosomes, segBAF.table, input_optimum_pair, dist_choice, distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA, gamma_param, read_depth, uninformative_BAF_threshold, allow100percent, reliabilityFile=NA, psi_min_initial=1.0, psi_max_initial=5.4, rho_min_initial=0.1, rho_max_initial=1.05, chr.names) # kjd 10-1-2014
{
  siglevel_BAF = 0.05 # kjd 21-2-2014
  # siglevel_BAF = 0.005 # kjd 21-2-2014
  
  maxdist_BAF = 0.01 # kjd 21-2-2014
  # maxdist_BAF = 0.005 # kjd 21-2-2014
  # maxdist_BAF = 0.001 # kjd 21-2-2014
  
  #siglevel_LogR = 0.05 # kjd 21-2-2014
  #maxdist_LogR = 0.1 # kjd 21-2-2014
  
  #DCW 160314 - much more lenient logR thresholds (allow anything!)
  siglevel_LogR = -0.01 # TODO: This parameter is pushed down to is.segment.clonal but not used there (maybe not used at all?)
  maxdist_LogR = 1 # TODO: This parameter is pushed down to is.segment.clonal but not used there (maybe not used at all?)

  
#   psi_min_initial = 1.0
#   psi_max_initial = 5.4
#   rho_min_initial = 0.1
#   rho_max_initial = 1.05
  
  ininitial_bounds = list( psi_min = psi_min_initial, psi_max = psi_max_initial, rho_min = rho_min_initial, rho_max = rho_max_initial )
  
  new_bounds = get_new_bounds( input_optimum_pair, ininitial_bounds ) # kjd 21-2-2014
  
  
  ch = chromosomes
  b = bafsegmented
  r = lrrsegmented[names(bafsegmented)]

  s = get_segment_info(lrrsegmented,segBAF.table)  
  # Make sure no segment of length 1 remains - TODO: this should not occur and needs to be prevented upstream
  s = s[s[,3] > 1,]
  dist_matrix_info <- create_distance_matrix_clonal( s, dist_choice, gamma_param, read_depth, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR, uninformative_BAF_threshold, new_bounds)# kjd 10-2-2013
  
  d = dist_matrix_info$distance_matrix # kjd 10-2-2013
  minimise = dist_matrix_info$minimise # kjd 10-2-2013
  
  #DCW 210314
  if(minimise){
	    best.distance = min(d)
  }else{
	  best.distance = max(d)
  }
  
  ref_seg_matrix = dist_matrix_info$ref_seg_matrix
  
  ref_major = dist_matrix_info$ref_major
  ref_minor = dist_matrix_info$ref_minor
  
  #########################################################
  
  ret = find_centroid_of_global_minima( d, ref_seg_matrix, ref_major, ref_minor, s, dist_choice, minimise, new_bounds, distancepng, gamma_param, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR, allow100percent, uninformative_BAF_threshold, read_depth) # kjd 28-2-2014
  optima_info_without_ref = ret$optima_info_without_ref
  optima_info = ret$optima_info
  
  nropt = optima_info$nropt
  psi_opt1 = optima_info$psi_opt1
  rho_opt1 = optima_info$rho_opt1
  ploidy_opt1 = optima_info$ploidy_opt1
  goodnessOfFit_opt1 = optima_info$goodnessOfFit_opt1

  distance.from.ref.seg = goodnessOfFit_opt1
  
  is.ref.better = F
  if (is.na(rho_opt1)) {
	  print("reference segment did not provide a possible solution")
  } else if(psi_opt1>= psi_min_initial & psi_opt1<= psi_max_initial & rho_opt1>= rho_min_initial & rho_opt1<= rho_max_initial & ((minimise & distance.from.ref.seg<best.distance)|(!minimise & distance.from.ref.seg>best.distance))){
	  is.ref.better = T
	  print("reference segment gives better results than grid search")
  } else {
	  print("reference segment gives no better results than grid search. Reverting to grid search solution")
  }

  psi_without_ref = optima_info_without_ref$psi_opt1
  rho_without_ref = optima_info_without_ref$rho_opt1
  ploidy_without_ref = optima_info_without_ref$ploidy_opt1
  goodnessOfFit_without_ref = optima_info_without_ref$goodnessOfFit_opt1
	
  #########################################################


  if(nropt>0) {
	
	#310314 DCW - always use grid search solution, because ref segment sometimes gives strange results
	#if(is.ref.better){
	#	rho = rho_opt1
	#	psi = psi_opt1
	#	ploidy = ploidy_opt1
	#	goodnessOfFit = goodnessOfFit_opt1
	#	print("ref segment gives best solution. Using this solution for plotting")
	#}else{
    rho = rho_without_ref
    psi = psi_without_ref
    ploidy = ploidy_without_ref
    goodnessOfFit = goodnessOfFit_without_ref*100
    #print("grid search gives best solution. Using this solution for plotting")
	#}
   
    nAfull = (rho-1-(b-1)*2^(r/gamma_param)*((1-rho)*2+rho*psi))/rho
    nBfull = (rho-1+b*2^(r/gamma_param)*((1-rho)*2+rho*psi))/rho
    nA = pmax(round(nAfull),0)
    nB = pmax(round(nBfull),0)
    
    rBacktransform = gamma_param*log((rho*(nA+nB)+(1-rho)*2)/((1-rho)*2+rho*psi),2)
    bBacktransform = (1-rho+rho*nB)/(2-2*rho+rho*(nA+nB))
    rConf = ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
    bConf = ifelse(bBacktransform!=0.5,pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))),NA)
    #DCW 150711 - get deviations from expected values
    if(!is.na(reliabilityFile)){
		  write.table(data.frame(segmentedBAF=b,backTransformedBAF=bBacktransform,confidenceBAF=bConf,segmentedR=r,backTransformedR=rBacktransform,confidenceR=rConf,nA=nA,nB=nB,nAfull=nAfull,nBfull=nBfull), reliabilityFile,sep=",",row.names=F)
	  }
    confidence = ifelse(is.na(rConf),bConf,ifelse(is.na(bConf),rConf,(rConf+bConf)/2))

	
  	# Make plots
  	if (!is.na(copynumberprofilespng)) { png(filename = copynumberprofilespng, width = 2000, height = 500, res = 200) }
  	ASCAT::ascat.plotAscatProfile(n1all = nA, n2all = nB, heteroprobes = TRUE, ploidy = ploidy, rho = rho, goodnessOfFit = goodnessOfFit, nonaberrant = FALSE, ch = ch, lrr = lrr, bafsegmented = bafsegmented, chrs=chr.names)
  	if (!is.na(copynumberprofilespng)) { dev.off() }
      
  	# separated plotting from logic: create nonrounded copy number profile plot here
  	if (!is.na(nonroundedprofilepng)) { png(filename = nonroundedprofilepng, width = 2000, height = 500, res = 200) }
  	ASCAT::ascat.plotNonRounded(ploidy = ploidy, rho = rho, goodnessOfFit = goodnessOfFit, nonaberrant = FALSE, nAfull = nAfull, nBfull = nBfull, bafsegmented = bafsegmented, ch = ch, lrr = lrr, chrs=chr.names)
  	if (!is.na(nonroundedprofilepng)) { dev.off() }
  }

  # Recalculate the psi_t for this rho using only clonal segments 
  psi_t = recalc_psi_t(psi_without_ref, rho_without_ref, gamma_param, lrrsegmented, segBAF.table, siglevel_BAF, maxdist_BAF, include_subcl_segments=F)

  # If there aren't any clonally fit segments, the above yields NA. In this case, revert to the original grid search psi_t
  if (is.na(psi_t)) {
	  print("Recalculated psi_t was NA, reverting to grid search solution. This occurs when no segment could be fit with a clonal state, check sample for contamination")
	  psi_t = psi_without_ref
  }
  
  output_optimum_pair = list(psi = psi_opt1, rho = rho_opt1, ploidy = ploidy_opt1)
  #output_optimum_pair_without_ref = list(psi = psi_without_ref, rho = rho_without_ref, ploidy = ploidy_without_ref)
  # Use the recalculated psi_t from the clonal segments as our final estimate of psi_t which is data driven with rho fixed
  output_optimum_pair_without_ref = list(psi = psi_t, rho = rho_without_ref, ploidy = ploidy_without_ref)
  return(list(output_optimum_pair=output_optimum_pair, output_optimum_pair_without_ref=output_optimum_pair_without_ref, distance = distance.from.ref.seg, distance_without_ref = best.distance, minimise = minimise, is.ref.better = is.ref.better)) # kjd 20-2-2014, adapted by DCW 140314
}

#' Recalculate psi_t based on rho and the available data
#' 
#' @param psi A psi estimate
#' @param rho A rho estimate
#' @param platform_gamma The platform specific LogR scaling parameter
#' @param lrrsegmented Segmented LogR, a vector with just the values
#' @param segBAF.table Segmented BAF, the full table
#' @param siglevel_BAF Significance level when testing wether a segment is clonal or subclonal given a rho/psi combination, parameter is used in \code{is.segment.clonal}
#' @param maxdist_BAF Max distance BAF is allowed to be away from the copy number solution before we don't trust the value and overrule a p-value, parameter required when determining the clonal status of a segment in \code{is.segment.clonal}
#' @param include_subcl_segments Boolean flag, supply TRUE if subclonal segments should be included when calculating psi_t, supply FALSE if only clonal segments should be included (default: TRUE)
#' @noRd
recalc_psi_t = function(psi, rho, gamma_param, lrrsegmented, segBAF.table, siglevel_BAF, maxdist_BAF, include_subcl_segments=T) {
  # Create segments of constant BAF/LogR  
  s = get_segment_info(lrrsegmented[rownames(segBAF.table)], segBAF.table)
  # Make sure no segment of length 1 remains - TODO: this should not occur and needs to be prevented upstream
  s = s[s[,3] > 1,]
  
  # Fetch all segments, if required check which ones are clonal with this rho/psi configuration
  segs = list()
  for (i in 1:nrow(s)) {
    read_depth = NA # Unused parameter
    maxdist_LogR = NA # Unused parameter
    siglevel_LogR = NA # Unused parameter
    segment_info = is.segment.clonal(LogR=s[i, "r"], 
                                     BAFreq=s[i, "b"], 
                                     BAF.length=s[i, "length"], 
                                     BAF.size=s[i, "size"], 
                                     BAF.mean=s[i, "mean"], 
                                     BAF.sd=s[i, "sd"], 
                                     read_depth=read_depth, 
                                     rho=rho, 
                                     psi=psi, 
                                     gamma_param=gamma_param, 
                                     siglevel_BAF=siglevel_BAF, 
                                     maxdist_BAF=maxdist_BAF, 
                                     siglevel_LogR=siglevel_LogR, 
                                     maxdist_LogR=maxdist_LogR)
    # Include this segment if we want to include all segments, or if we don't want subclonal segments include it only if its clonal
    if (include_subcl_segments | segment_info$is.clonal) {
      nMaj = segment_info$nMaj.test
      nMin = segment_info$nMin.test
      psi_t = calc_psi_t(nMaj+nMin, s[i, "r"], rho, gamma_param)
      segs[[length(segs)+1]] = data.frame(nMaj=nMaj, nMin=nMin, length=s[i, "length"], psi_t=psi_t)
    }
  }
  segs = do.call(rbind, segs)
      
  # Calculate psi_t as the weighted average copy number across all segments
  psi_t = sum(segs$psi_t * segs$length, na.rm=T) / sum(segs$length, na.rm=T)
  return(psi_t)
}


#' Calculate psi based on a reference segment and its associated logr
#' 
#' @param total_cn Integer representing the total clonal copynumber (i.e. nMajor+nMinor)
#' @param r The LogR of the segment with the total_cn copy number
#' @param rho A cellularity estimate
#' @param gamma_param Platform gamma parameter
#' @author sd11
#' @export
calc_psi_t = function(total_cn, r, rho, gamma_param) {
  psi = (rho*(total_cn)+2-2*rho)/(2^(r/gamma_param))
  psi_t = (psi-2*(1-rho))/rho
  return(psi_t)
}
