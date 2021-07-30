#PCF-ALGORITHM (KL):
### EXACT version
exactPcf <- function(y, kmin=5, gamma, yest) {
## Implementaion of exact PCF by Potts-filtering
	## x: input array of (log2) copy numbers
	## kmin: Mininal length of plateaus
	## gamma: penalty for each discontinuity
  	N <- length(y)
  	yhat <- rep(0,N);
  	if (N < 2*kmin) {
	     if (yest) {
 		     return(list(Lengde = N, sta = 1, mean = mean(y), nIntervals=1, yhat=rep(mean(y),N)))
       } else {
 		     return(list(Lengde = N, sta = 1, mean = mean(y), nIntervals=1))
 	     }
  	}
  	initSum <- sum(y[1:kmin])
  	initKvad <- sum(y[1:kmin]^2)
  	initAve <- initSum/kmin;
  	bestCost <- rep(0,N)
  	bestCost[kmin] <- initKvad - initSum*initAve
  	bestSplit <- rep(0,N)
  	bestAver <- rep(0,N)
  	bestAver[kmin] <- initAve
  	Sum <- rep(0,N)
  	Kvad <- rep(0,N)
  	Aver <- rep(0,N)
  	Cost <- rep(0,N)
  	kminP1=kmin+1
  	for (k in (kminP1):(2*kmin-1)) {
    		Sum[kminP1:k]<-Sum[kminP1:k]+y[k]
    		Aver[kminP1:k] <- Sum[kminP1:k]/((k-kmin):1)
    		Kvad[kminP1:k] <- Kvad[kminP1:k]+y[k]^2
    		bestAver[k] <- (initSum+Sum[kminP1])/k
    		bestCost[k] <- (initKvad+Kvad[kminP1])-k*bestAver[k]^2
  	}
  	for (n in (2*kmin):N) {
   		yn <- y[n]
   		yn2 <- yn^2
   		Sum[kminP1:n] <- Sum[kminP1:n]+yn
   		Aver[kminP1:n] <- Sum[kminP1:n]/((n-kmin):1)
   		Kvad[kminP1:n] <- Kvad[kminP1:n]+yn2
   		nMkminP1=n-kmin+1
   		Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)]+Kvad[kminP1:nMkminP1]-Sum[kminP1:nMkminP1]*Aver[kminP1:nMkminP1]+gamma
   		Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
   		cost <- Cost[Pos]
   		aver <- Aver[Pos]
   		totAver <- (Sum[kminP1]+initSum)/n
   		totCost <- (Kvad[kminP1]+initKvad) - n*totAver*totAver
   		if (totCost < cost) {
       		Pos <- 1
          cost <- totCost
          aver <- totAver
   		}
   		bestCost[n] <- cost
   		bestAver[n] <- aver
   		bestSplit[n] <- Pos-1
 	}
 	n <- N
	antInt <- 0
	if(yest){
 		while (n > 0) {
   			yhat[(bestSplit[n]+1):n] <- bestAver[n]
   			n <- bestSplit[n]
			antInt <- antInt+1
 		}
 	} else {
 		while (n > 0) {
	   		n <- bestSplit[n]
			antInt <- antInt+1
 		}
 	}
	n <- N  #nProbes   Spr Knut, fant ikke nProbes noe sted..
	lengde <- rep(0,antInt)
	start <- rep(0,antInt)
	verdi <- rep(0,antInt)
	oldSplit  <- n
	antall <- antInt
	while (n > 0) {
    start[antall] <- bestSplit[n]+1
		lengde[antall] <- oldSplit-bestSplit[n]
		verdi[antall] <- bestAver[n]
		n <- bestSplit[n]
		oldSplit <- n
		antall <- antall-1
 	}
	if (yest) {
 		return(list(Lengde = lengde, sta = start, mean = verdi, nIntervals=antInt, yhat=yhat))
 	} else {
 		return(list(Lengde = lengde, sta = start, mean = verdi, nIntervals=antInt))
 	}
}



selectFastPcf <- function(x,kmin,gamma,yest){
	xLength <- length(x)
	if (xLength< 1000) {
		result<-runFastPcf(x,kmin,gamma,0.15,0.15,yest)
	} else {
	if (xLength < 15000){
		result<-runFastPcf(x,kmin,gamma,0.12,0.05,yest)
	} else  {
		result<-runPcfSubset(x,kmin,gamma,0.12,0.05,yest)
		}
	}
	return(result)
}


runFastPcf <- function(x,kmin,gamma,frac1,frac2,yest){
	antGen <- length(x)
	mark<-filterMarkS4(x,kmin,8,1,frac1,frac2,0.02,0.9)
	mark[antGen]=TRUE
	dense <- compact(x,mark)
	#print(dense$Nr)
	#print(frac2)
	result<-PottsCompact(kmin,gamma,dense$Nr,dense$Sum,dense$Sq,yest)
	return(result)
}

runPcfSubset <- function(x,kmin,gamma,frac1,frac2,yest){
	SUBSIZE <- 5000
	antGen <- length(x)
	mark<-filterMarkS4(x,kmin,8,1,frac1,frac2,0.02,0.9)
	markInit<-c(mark[1:(SUBSIZE-1)],TRUE)
	compX<-compact(x[1:SUBSIZE],markInit)
	mark2 <- rep(FALSE,antGen)
	mark2[1:SUBSIZE] <- markWithPotts(kmin,gamma,compX$Nr,compX$Sum,compX$Sq,SUBSIZE)
	mark2[4*SUBSIZE/5]<-TRUE
	start <- 4*SUBSIZE/5+1
	while(start + SUBSIZE < antGen){
		slutt<-start+SUBSIZE-1
		markSub<-c(mark2[1:(start-1)],mark[start:slutt])
		markSub[slutt] <- TRUE
		compX<-compact(x[1:slutt],markSub)
		mark2[1:slutt] <- markWithPotts(kmin,gamma,compX$Nr,compX$Sum,compX$Sq,slutt)
		start <- start+4*SUBSIZE/5
		mark2[start-1]<-TRUE
	}
	markSub<-c(mark2[1:(start-1)],mark[start:antGen])
	compX<-compact(x,markSub)
	result <- PottsCompact(kmin,gamma,compX$Nr,compX$Sum,compX$Sq,yest)
      	return(result)
}

PottsCompact <- function(kmin, gamma, nr, res, sq, yest) {
## Potts filtering on compact array;
	## kmin: minimal length of plateau
	## gamma: penalty for discontinuity
	## nr: number of values between breakpoints
	## res: sum of values between breakpoints
	## sq: sum of squares of values between breakpoints

  	N <- length(nr)
  	Ant <- rep(0,N)
  	Sum <- rep(0,N)
  	Kvad <- rep(0,N)
  	Cost <- rep(0,N)
  	if (sum(nr) < 2*kmin){
  		estim <- sum(res)/sum(nr)
		return(estim)
  	}
  	initAnt <- nr[1]
  	initSum <- res[1]
  	initKvad <- sq[1]
  	initAve <- initSum/initAnt
  	bestCost <- rep(0,N)
  	bestCost[1] <- initKvad - initSum*initAve
  	bestSplit <- rep(0,N)
  	k <- 2
  	while(sum(nr[1:k]) < 2*kmin) {
		Ant[2:k] <- Ant[2:k]+nr[k]
    		Sum[2:k]<-Sum[2:k]+res[k]
    		Kvad[2:k] <- Kvad[2:k]+sq[k]
    		bestCost[k] <- (initKvad+Kvad[2])-(initSum+Sum[2])^2/(initAnt+Ant[2])
		k <- k+1	
  	}
  	for (n in k:N) {
		Ant[2:n] <- Ant[2:n]+nr[n]
   		Sum[2:n] <- Sum[2:n]+res[n]
   		Kvad[2:n] <- Kvad[2:n]+sq[n]
		limit <- n
		while(limit > 2 & Ant[limit] < kmin) {limit <- limit-1}
   		Cost[2:limit] <- bestCost[1:limit-1]+Kvad[2:limit]-Sum[2:limit]^2/Ant[2:limit]
   		Pos <- which.min(Cost[2:limit])+ 1
   		cost <- Cost[Pos]+gamma
   		totCost <- (Kvad[2]+initKvad) - (Sum[2]+initSum)^2/(Ant[2]+initAnt)
   		if (totCost < cost) {
       			Pos <- 1
       			cost <- totCost
   		}
   		bestCost[n] <- cost
   		bestSplit[n] <- Pos-1
 	}
  	if (yest) {
		yhat<-rep(0,N)
		res<-findEst(bestSplit,N,nr,res,TRUE)
	} else {
		res<-findEst(bestSplit,N,nr,res,FALSE)
	}
	return(res)
}

compact <- function(y,mark){
  ## accumulates numbers of observations, sums and 
  ## sums of squares between potential breakpoints
	return(list(
	        Nr = diff(append(0, which(mark))),
	        Sum = diff(append(0, cumsum(y)[mark])),
	        Sq = diff(append(0, cumsum(y ^ 2)[mark]))))
}

findEst <- function(bestSplit,N,Nr,Sum,yest){
	n<-N
	lengde<-rep(0,N)
	antInt<-0
	while (n>0){
		antInt<-antInt+1
		lengde[antInt] <- n-bestSplit[n]
		n<-bestSplit[n]
	}
	lengde<-lengde[antInt:1]
	lengdeOrig<-rep(0,antInt)
	startOrig<-rep(1,antInt+1)
	verdi<-rep(0,antInt)
	start<-rep(1,antInt+1)
	for(i in 1:antInt){
		start[i+1] <- start[i]+lengde[i]
		lengdeOrig[i] <- sum(Nr[start[i]:(start[i+1]-1)])
		startOrig[i+1] <- startOrig[i]+lengdeOrig[i]
		verdi[i]<-sum(Sum[start[i]:(start[i+1]-1)])/lengdeOrig[i]
	}
	
	if(yest){
		yhat<-rep(0,startOrig[antInt+1]-1)
		for (i in 1:antInt){
			yhat[startOrig[i]:(startOrig[i+1]-1)]<-verdi[i]
		}
		startOrig<-startOrig[1:antInt]
		return(list(Lengde=lengdeOrig,sta=startOrig,mean=verdi,nIntervals=antInt,yhat=yhat))
	} else {
		startOrig<-startOrig[1:antInt]
		return(list(Lengde=lengdeOrig,sta=startOrig,mean=verdi,nIntervals=antInt))
	}
	
}


markWithPotts <- function(kmin, gamma, nr, res, sq, subsize) {
## Potts filtering on compact array;
	## kmin: minimal length of plateau
	## gamma: penalty for discontinuity
	## nr: number of values between breakpoints
	## res: sum of values between breakpoints
	## sq: sum of squares of values between breakpoints

  	N <- length(nr)
  	Ant <- rep(0,N)
  	Sum <- rep(0,N)
  	Kvad <- rep(0,N)
  	Cost <- rep(0,N)
	markSub <- rep(FALSE,N)
  	initAnt <- nr[1]
  	initSum <- res[1]
  	initKvad <- sq[1]
  	initAve <- initSum/initAnt
  	bestCost <- rep(0,N)
  	bestCost[1] <- initKvad - initSum*initAve
  	bestSplit <- rep(0,N)
  	k <- 2
  	while(sum(nr[1:k]) < 2*kmin) {
		Ant[2:k] <- Ant[2:k]+nr[k]
    		Sum[2:k]<-Sum[2:k]+res[k]
    		Kvad[2:k] <- Kvad[2:k]+sq[k]
    		bestCost[k] <- (initKvad+Kvad[2])-(initSum+Sum[2])^2/(initAnt+Ant[2])
		k <- k+1	
  	}
  	for (n in k:N) {
		Ant[2:n] <- Ant[2:n]+nr[n]
   		Sum[2:n] <- Sum[2:n]+res[n]
   		Kvad[2:n] <- Kvad[2:n]+sq[n]
		limit <- n
		while(limit > 2 & Ant[limit] < kmin) {limit <- limit-1}
   		Cost[2:limit] <- bestCost[1:limit-1]+Kvad[2:limit]-Sum[2:limit]^2/Ant[2:limit]
   		Pos <- which.min(Cost[2:limit])+ 1
   		cost <- Cost[Pos]+gamma
   		totCost <- (Kvad[2]+initKvad) - (Sum[2]+initSum)^2/(Ant[2]+initAnt)
   		if (totCost < cost) {
       			Pos <- 1
       			cost <- totCost
   		}
   		bestCost[n] <- cost
   		bestSplit[n] <- Pos-1
		markSub[Pos-1] <- TRUE
 	}
 	help<-findMarks(markSub,nr,subsize)
	return(help=help)
}


findMarks <- function(markSub,Nr,subsize){
	## markSub: marks in compressed scale
	## NR: number of observations between potenstial breakpoints
	mark<-rep(FALSE,subsize)  ## marks in original scale
	if(sum(markSub)<1) {return(mark)} else {	
		N<-length(markSub)
		ant <- seq(1:N)
		help <- ant[markSub]
		lengdeHelp<-length(help)
		help0 <- c(0,help[1:(lengdeHelp-1)])
		lengde <- help-help0
		start<-1
		oldStart<-1
		startOrig<-1
		for(i in 1:lengdeHelp){
			start <- start+lengde[i]
			lengdeOrig <- sum(Nr[oldStart:(start-1)])
			startOrig <- startOrig+lengdeOrig
			mark[startOrig-1]<-TRUE
			oldStart<-start
		}
		return(mark)
	}
	
}


compact <- function(y,mark){
## accumulates numbers of observations, sums and 
## sums of squares between potential breakpoints
## y:  array to be compacted
## mark:  logical array of potential breakpoints
	tell<-seq(1:length(y))
	cCTell<-tell[mark]
	Ncomp<-length(cCTell)
	lowTell<-c(0,cCTell[1:(Ncomp-1)])
	ant<-cCTell-lowTell
	cy<-cumsum(y)
	cCcy<-cy[mark]
	lowcy<-c(0,cCcy[1:(Ncomp-1)])
	sum<-cCcy-lowcy
	cy2<-cumsum(y^2)
	cCcy2<-cy2[mark]
	lowcy2<-c(0,cCcy2[1:(Ncomp-1)])
	sq<-cCcy2-lowcy2
	return(list(Nr=ant,Sum=sum,Sq=sq))
}

filterMarkS4 <- function(x,kmin,L,L2,frac1,frac2,frac3,thres){
## marks potential breakpoints, partially by a two 6*L and 6*L2 highpass
## filters (L>L2), then by a filter seaching for potential kmin long segments
    lengdeArr <- length(x)
    xc<-cumsum(x)
	xc<-c(0,xc)
	ind11<-1:(lengdeArr-6*L+1)
	ind12<-ind11+L
	ind13<-ind11+3*L
	ind14<-ind11+5*L
	ind15<-ind11+6*L
	cost1<-abs(4*xc[ind13]-xc[ind11]-xc[ind12]-xc[ind14]-xc[ind15])	
	cost1<-c(rep(0,3*L-1),cost1,rep(0,3*L))
	##mark shortening in here
	in1<-1:(lengdeArr-6)
	in2<-in1+1
	in3<-in1+2
	in4<-in1+3
	in5<-in1+4
	in6<-in1+5
	in7<-in1+6
	test<-pmax(cost1[in1],cost1[in2],cost1[in3],cost1[in4],cost1[in5],cost1[in6],cost1[in7])
	test<-c(rep(0,3),test,rep(0,3))
	cost1B<-cost1[cost1>=thres*test]
	frac1B<-min(0.8,frac1*length(cost1)/length(cost1B))
	limit <- quantile(cost1B,(1-frac1B),names=FALSE)
	mark<-(cost1>limit)&(cost1>0.9*test)	
	
	
	ind21<-1:(lengdeArr-6*L2+1)
	ind22<-ind21+L2
	ind23<-ind21+3*L2
	ind24<-ind21+5*L2
	ind25<-ind21+6*L2
	cost2<-abs(4*xc[ind23]-xc[ind21]-xc[ind22]-xc[ind24]-xc[ind25])
	limit2 <- quantile(cost2,(1-frac2),names=FALSE)
	mark2<-(cost2>limit2)
	mark2<-c(rep(0,3*L2-1),mark2,rep(0,3*L2))
	if(3*L>kmin){
		mark[kmin:(3*L-1)]<-TRUE
		mark[(lengdeArr-3*L+1):(lengdeArr-kmin)]<-TRUE
	}
	else
	{
		mark[kmin]<- TRUE
		mark[lengdeArr-kmin]<-TRUE
	}

	if(kmin>1){
		ind1<-1:(lengdeArr-3*kmin+1)
		ind2<-ind1+3*kmin
		ind3<-ind1+kmin
		ind4<-ind1+2*kmin
		shortAb <- abs(3*(xc[ind4]-xc[ind3])-(xc[ind2]-xc[ind1]))
		in1<-1:(length(shortAb)-6)
		in2<-in1+1
		in3<-in1+2
		in4<-in1+3
		in5<-in1+4
		in6<-in1+5
		in7<-in1+6
		test<-pmax(shortAb[in1],shortAb[in2],shortAb[in3],shortAb[in4],shortAb[in5],shortAb[in6],shortAb[in7])
		test<-c(rep(0,3),test,rep(0,3))
		cost1C<-shortAb[shortAb>=thres*test]
		frac1C<-min(0.8,frac3*length(shortAb)/length(cost1C))
		limit3 <- quantile(cost1C,(1-frac1C),names=FALSE)
		markH1<-(shortAb>limit3)&(shortAb>thres*test)
		markH2<-c(rep(FALSE,(kmin-1)),markH1,rep(FALSE,2*kmin))
		markH3<-c(rep(FALSE,(2*kmin-1)),markH1,rep(FALSE,kmin))
		mark<-mark|mark2|markH2|markH3
	} else {
		mark<-mark|mark2
	}

	if(3*L>kmin){
		mark[1:(kmin-1)]<-FALSE
		mark[kmin:(3*L-1)]<-TRUE
		mark[(lengdeArr-3*L+1):(lengdeArr-kmin)]<-TRUE
		mark[(lengdeArr-kmin+1):(lengdeArr-1)]<-FALSE
		mark[lengdeArr]<-TRUE	
	}
	else
	{
		mark[1:(kmin-1)]<-FALSE
		mark[(lengdeArr-kmin+1):(lengdeArr-1)]<-FALSE
		mark[lengdeArr]<-TRUE
		mark[kmin]<- TRUE
		mark[lengdeArr-kmin]<-TRUE
	}

	return(mark)
}

#Get mad SD-estimate

##Input:
### x: vector of observations for which mad Sd is to be calculated
### k: window size to be used in median filtering

##Output:
### SD: mad sd estimate

##Required by:
### multiPcf
### fastPcf
### pcf
### aspcf


##Requires:
### medianFilter




getMad <- function(x,k=25){
  
  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
  x <- x[x!=0]
  
  #Calculate runMedian  
  runMedian <- medianFilter(x,k)
  
  dif <- x-runMedian
  SD <- mad(dif)
 
	return(SD)
}


#########################################################################
# Function to calculate running median for a given a window size
#########################################################################

##Input:
### x: vector of numeric values
### k: window size to be used for the sliding window (actually half-window size)

## Output:
### runMedian : the running median corresponding to each observation

##Required by:
### getMad
### medianFilter


##Requires:
### none

medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1
  
  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }
  
  runMedian <- runmed(x,k=filtWidth,endrule="median")

  return(runMedian)

}
