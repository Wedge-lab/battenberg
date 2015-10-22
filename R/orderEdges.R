#' Convenience function that orders edges or squares
#' @author dw9, kd7
#' @noRd
orderEdges = function(levels, l, ntot,x,y) {
  nMaj1 = NULL
  nMin1 = NULL
  nMaj2 = NULL
  nMin2 = NULL
  
  # case 1 or 2a:
  if(l>levels[3]) {
    #LogR criterion: ntot < x+y+1
    if(ntot < x+y+1) {
      # take the six options, sorted according to LogR priority (3+3) + simplicity (1+2+1+2)
      nMaj1 = c(y,y-1,y,
                y+1,y+1,y+1)
      nMin1 = c(x,x,x,
                x,x-1,x)
      nMaj2 = c(y+1,y+1,y+2,
                y+1,y+1,y+1)
      nMin2 = c(x,x,x,
                x+1,x+1,x+2)
    }
    else {
      nMaj1 = c(y+1,y+1,y+1,
                y,y-1,y)
      nMin1 = c(x,x-1,x,
                x,x,x)
      nMaj2 = c(y+1,y+1,y+1,
                y+1,y+1,y+2)
      nMin2 = c(x+1,x+1,x+2,
                x,x,x)
    }
  }
  # case 2c:
  else if (l>levels[2]) {
    if(ntot < x+y+1) {
      nMaj1 = c(y,y,y,
                y+1,y+1,y+1)
      nMin1 = c(x,x-1,x,
                x,x-1,x)
      nMaj2 = c(y,y,y,
                y+1,y+1,y+1)
      nMin2 = c(x+1,x+1,x+2,
                x+1,x+1,x+2)
    }
    else {
      nMaj1 = c(y+1,y+1,y+1,
                y,y,y)
      nMin1 = c(x,x-1,x,
                x,x-1,x)
      nMaj2 = c(y+1,y+1,y+1,
                y,y,y)
      nMin2 = c(x+1,x+1,x+2,
                x+1,x+1,x+2)
    }
  }
  # case 2b:
  else {
    if(ntot < x+y+1) {
      nMaj1 = c(y,y,y,
                y,y-1,y)
      nMin1 = c(x,x-1,x,
                x+1,x+1,x+1)
      nMaj2 = c(y,y,y,
                y+1,y+1,y+2)
      nMin2 = c(x+1,x+1,x+2,
                x+1,x+1,x+1)
    }
    else {
      nMaj1 = c(y,y-1,y,
                y,y,y)
      nMin1 = c(x+1,x+1,x+1,
                x,x-1,x)
      nMaj2 = c(y+1,y+1,y+2,
                y,y,y)
      nMin2 = c(x+1,x+1,x+1,
                x+1,x+1,x+2)
    }
  }
  #DCW 260314 - avoid negative CNs
  negative.CN = which(nMaj1<0|nMin1<0|nMaj2<0|nMin2<0)
  if(length(negative.CN)>0){
    nMaj1[negative.CN]=NA
    nMin1[negative.CN]=NA
    nMaj2[negative.CN]=NA
    nMin2[negative.CN]=NA    	
    return(cbind(nMaj1,nMin1,nMaj2,nMin2))
  }else{
    return(cbind(nMaj1,nMin1,nMaj2,nMin2))
  }
}


#' Function that fetches the nearest edge for a given a rho, psi, BAF and major and minor allele
#' that corresponds to a certain mix of two copy number states. It first identifies the nearest edge 
#' and then just compares the vertices at the end of this edge to find the best corner.
#' @author dw9, kd7
#' @noRd
GetNearestCorners_bestOption <-function( rho, psi, BAFreq, nMajor, nMinor ) {	
  nMaj = c(floor(nMajor),ceiling(nMajor),floor(nMajor),ceiling(nMajor))
  nMin = c(ceiling(nMinor),ceiling(nMinor),floor(nMinor),floor(nMinor))
  x = floor(nMinor)
  y = floor(nMajor)
  
  # total copy number, to determine priority options
  ntot = nMajor + nMinor
  
  BAF_levels = (1-rho+rho*nMaj)/(2-2*rho+rho*(nMaj+nMin))
  #problem if rho=1 and nMaj=0 and nMin=0
  BAF_levels[nMaj==0 & nMin==0] = 0.5
  
  nMaj1 = NULL
  nMin1 = NULL
  nMaj2 = NULL
  nMin2 = NULL
  
  # case 1 or 2a:
  #if( is.finite(BAF_levels[3]) && (BAFreq>BAF_levels[3]) ) { # kjd 14-2-2014
  if(BAFreq>BAF_levels[3]) { #DCW
    #LogR criterion: ntot < x+y+1
    if(ntot < x+y+1) {
      # take the six options, sorted according to LogR priority (3+3) + simplicity (1+2+1+2)
      nMaj1 = y
      nMin1 = x
      nMaj2 = y+1
      nMin2 = x
    }
    else {
      nMaj1 = y+1
      nMin1 = x
      nMaj2 = y+1
      nMin2 = x+1
    }
  }
  # case 2c:
  #else if( is.finite(BAF_levels[2]) &&  (BAFreq>BAF_levels[2]) ) { # kjd 14-2-2014
  else if(BAFreq>BAF_levels[2]) { #DCW
    if(ntot < x+y+1) {
      nMaj1 = y
      nMin1 = x
      nMaj2 = y
      nMin2 = x+1
    }
    else {
      nMaj1 = y+1
      nMin1 = x
      nMaj2 = y+1
      nMin2 = x+1
    }
  }
  # case 2b:
  else {
    if(ntot < x+y+1) {
      nMaj1 = y
      nMin1 = x
      nMaj2 = y
      nMin2 = x+1
    }
    else {
      nMaj1 = y
      nMin1 = x+1
      nMaj2 = y+1
      nMin2 = x+1
    }
  }
  
  nMaj_vect = c( nMaj1, nMaj2 )
  nMin_vect = c( nMin1, nMin2 )
  
  nearest_segment = list( nMaj = nMaj_vect, nMin = nMin_vect )
  
  return( nearest_segment )
}
