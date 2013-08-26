##################################################
#			segmentation functions
##################################################
simSegData<-function(nseg=10, J=3, soj, emis, seed=1234, toPlot=FALSE){
	# focus on the soj settings,
	# keep the difference between emission of different states small

	# pool of segments for each states to sample from
	#check J to see if it conforms with soj and emis
	if(!is.null(soj)){
		if(! soj$type %in% c('gamma', 'pois', 'nbinom')){
			stop("soj type not supported")
		}
		paraLen<-sapply(soj, length)
		if(! all(paraLen[names(paraLen)!='type']==J)){
			stop("incorrect length for the soj parameter")
		}
	}
	if(!is.null(emis)){
		if(! emis$type %in% c('norm', 't', 'pois','nbinom')){
			stop("emis type not supported")
		}
		paraLen<-sapply(emis, length)
		if(!all(paraLen[names(paraLen)!='type']==J)){
			stop("incorrect length for the emis parameter")
		}
	}	

	# generate pool of segments, with soj related length, and states related signal
	segLen<-sapply(1:J, function(j) simUvSegData(nseg, j, soj, seed))
	segExp<-sapply(1:J, function(j) sapply(segLen[,j], function(l) simUvSegData(l, j, emis)))
	
	sel<-sample(1:J, size=nseg, replace=T)
	states<-runValue(Rle(sel))
	nsel<-length(states)
	
	L<-segLen[(states-1)*nseg+1:nsel]
	E<-unlist(segExp[(states-1)*nseg+1:nsel])
	
	res<-list(L=L,E=E,S=states)
	if(toPlot){
		filename<-paste('Simulation.', gsub(':', '-', gsub(' ', '_', date())) ,'.pdf', sep='')
		pdf(filename)
		ts.plot(E)
		dev.off()
		res<-c(res, pdf=filename)
	}
	return(res)
}


simUvSegData<-function(n, j, param, seed=NULL){
	if(!is.null(seed))	set.seed(seed)
	x<-switch(param$type,
		norm = rnorm(n, mean=param$mean[j], sd=param$sd[j]),
		t = rt(n, ncp=param$ncp[j], df=param$df[j]),
		gamma = rgamma(n, shape=param$shape[j], scale=param$scale[j]),
		pois = rpois(n, lambda=param$lambda[j]),
		nbinom = rnbinom(n, mu=param$mu[j], size=param$size[j]))
}


mat2list<-function(x, nr=NULL){
	if(!is.null(nr)) x<-matrix(x, nrow=nr)
	L = vector(mode="list", length=nrow(x))
	for(i in 1:nrow(x)) {
		L[[i]] = cbind(estimate=x[i, 0:(i-1)]) 
	}
	return(L)
}


preClustGrp<-function(x, grp=NULL, cluster.m=NULL, deepSplit=4, minClusterSize=4){
	clustmethods<-c('ward','single','complete','average','mcquitty','median','centroid')
	nc<-ncol(x)
	if (!is.null(grp)) {
		if(length(grp)!=nc) 
			stop('length of grp differs from number of column in x !!!')
		if(!is.null(cluster.m)) 
			warning('clustering would not be performed, using the grp input instead !!!')	
		grp<-as.character(grp)
	}	else if(!is.null(cluster.m) && nc >= 4){
		# no grp and has cluster.m
		if(length(find.package('dynamicTreeCut', quiet=T))==0) {
			warning("'dynamicTreeCut' is not found, clustering result ignored, treated as one group!!!")
			grp<-rep('1', nc)
		} else{
			cluster.m<-match.arg(cluster.m, clustmethods)
			DM <- dist(t(x))
			HC <- hclust(DM, method=cluster.m)
			require(dynamicTreeCut)
			grp<-cutreeDynamic(dendro = HC, distM = as.matrix(DM), deepSplit = deepSplit, pamRespectsDendro = FALSE, minClusterSize = minClusterSize, verbose=0) 
		}
	} else {
		grp<-rep('1', nc)
	}
	return(grp)
}

## main function 
biomvRmgmr<-function(x, xPos=NULL, xRange=NULL, usePos='start', cutoff=NULL, q=0.9, high=TRUE, minrun=5, maxgap=2, splitLen=Inf, poolGrp=FALSE, grp=NULL, cluster.m=NULL, avg.m='median', trim=0, na.rm=TRUE){
	
	if (!is.numeric(x) &&  !is.matrix(x) && class(x)!='GRanges') 
        stop("'x' must be a numeric vector or matrix or a GRanges object.")
    if(class(x)=='GRanges') {
    	xid<-names(values(x))
    	xRange<-x
    	mcols(xRange)<-NULL
    	x<-as.matrix(values(x))
		if( any(sapply(unique(seqnames(xRange)), function(s) length(unique(strand(xRange[seqnames(xRange)==s]))))!=1))
    		stop('For some sequence, there are data appear on both strands !')
    } else if(length(dim(x))==2){
		xid<-colnames(x)
		x<-as.matrix(x)
	} else {
		message('No dim attributes, coercing x to a matrix with 1 column.')
		x <- matrix(as.numeric(x), ncol=1)
	}
	nr<-nrow(x) 
	nc<-ncol(x)
	if(is.null(xid)){
		xid<-paste('S', seq_len(nc), sep='')
		colnames(x)<-xid
	}
	
	## some checking on xpos and xrange, xrange exist then xpos derived from xrange,
	if(!is.null(xRange) && (class(xRange)=='GRanges' || class(xRange)=='IRanges') && !is.null(usePos) && length(xRange)==nr && usePos %in% c('start', 'end', 'mid')){
		if(usePos=='start'){
			xPos<-start(xRange)
		} else if(usePos=='end'){
			xPos<-end(xRange)
		} else {
			xPos<-(start(xRange)+end(xRange))/2
		}
	} else {
		# no valid xRange, set it to null
		message('no valid xRange and usePos found, check if you have specified xRange / usePos.')
		xRange<- NULL
	} 
	if (is.null(xPos) || !is.numeric(xPos) || length(xPos)!=nr){
		message("No valid positional information found. Re-check if you have specified any xPos / xRange.")
		xPos<-NULL
	}
 
    # check grp setting, cluster if needed, otherwise treat as one group	
   	if(!is.null(grp)) grp<-as.character(grp)
	grp<-preClustGrp(x, grp=grp, cluster.m=cluster.m)
	
	## build xRange if not a GRanges for the returning object
	if(is.null(xRange) || class(xRange) != 'GRanges'){
		if(!is.null(xRange) && class(xRange) == 'IRanges'){
			xRange<-GRanges(seqnames='sampleseq', xRange)	
		} else 	if(!is.null(xPos)){
			xRange<-GRanges(seqnames='sampleseq', IRanges(start=xPos, width=1))	
		} else {
			xRange<-GRanges(seqnames='sampleseq', IRanges(start=seq_len(nr), width=1))	
		}
	}
	# get seqnames status	
	seqs<-unique(as.character(seqnames(xRange)))
	## initialize the output vectors
	res<-GRanges()
	seqlevels(res)<-seqlevels(xRange)
	# we have more than one seq to batch
	for(s in seq_along(seqs)){
		message(sprintf("Processing sequence %s\n", seqs[s]))
		r<-which(as.character(seqnames(xRange)) == seqs[s])
		
		for(g in unique(grp)){
			gi<-grp==g
			message(sprintf("Building segmentation model for group %s ...", g))
			# fixme, worth thinking, if columns within one group should be pooled.
			if(poolGrp){
				Ilist<-maxGapminRun(x=apply(as.matrix(x[r, gi]), 1, median, na.rm=na.rm), xRange=ranges(xRange)[r], cutoff=cutoff, q=q, high=high, minrun=minrun, maxgap=maxgap, splitLen=splitLen, na.rm=na.rm)
				if(length(Ilist$IS)>0){
					tores<-GRanges(seqnames=as.character(seqs[s]), 
								IRanges(start=rep(start(xRange)[r][Ilist$IS], sum(gi)), end=rep(end(xRange)[r][Ilist$IE], sum(gi))), 
								strand=strand(xRange)[r][Ilist$IS],
								SAMPLE=rep(xid[gi], each=length(Ilist$IS)), 
								STATE=rep('HI', length(Ilist$IS)*sum(gi)), 
								AVG=sapply(which(gi), function(c) sapply(seq_along(Ilist$IS), function(t) avgFunc(x[Ilist$IS[t]:Ilist$IE[t],c], avg.m=avg.m, trim=trim, na.rm=na.rm)))
					)
					seqlevels(tores)<-seqlevels(xRange)
					res<-c(res, tores)					
				}
			} else {
				for(c in which(gi)){
				Ilist<-maxGapminRun(x=x[r,c], xRange=ranges(xRange)[r], cutoff=cutoff, q=q, high=high, minrun=minrun, maxgap=maxgap, splitLen=splitLen, na.rm=na.rm)
				if(length(Ilist$IS)>0){
					tores<-GRanges(seqnames=as.character(seqs[s]), 
						IRanges(start=start(xRange)[r][Ilist$IS], end=end(xRange)[r][Ilist$IE]), 
						strand=strand(xRange)[r][Ilist$IS],
						SAMPLE=rep(xid[c], length(Ilist$IS)), 
						STATE=rep(ifelse(high, 'HIGH', 'LOW'), length(Ilist$IS)), 
						AVG=as.numeric(sapply(seq_along(Ilist$IS),  function(t) avgFunc(x[Ilist$IS[t]:Ilist$IE[t],c], avg.m=avg.m, trim=trim, na.rm=na.rm)))
					)
					seqlevels(tores)<-seqlevels(xRange)
					res<-c(res, tores)	
				}
			}
			
			} # end c for
			message(sprintf("Building segmentation model for group %s complete.", g))
		} # end for g
		message(sprintf("Processing sequence %s complete.", seqs[s]))
	} # end for s

	values(xRange)<-DataFrame(x,  row.names = NULL)
	new("biomvRCNS",  
		x = xRange, res = res,
		param=list(maxgap=maxgap, minrun=minrun, q=q, cutoff=cutoff, splitLen=splitLen, grp=grp, cluster.m=cluster.m, poolGrp=poolGrp, avg.m=avg.m, trim=trim, na.rm=na.rm)
	)
}


##################################################
#  maxgap minrun algo, utilizing RLE, mostly for transcript detection...
##################################################
maxGapminRun<-function(x, xPos=NULL, xRange=NULL, cutoff=NULL, q=0.9, high=TRUE, minrun=5, maxgap=2, splitLen=Inf, na.rm=TRUE){
	# optional position vector, if not present, use index as postion
	# output values, index of active intervals, plus func parameters 
	
	if(is.null(cutoff)){
		cutoff<- as.numeric(quantile(x, q, na.rm=na.rm))
	} else if (cutoff >= max(x, na.rm=na.rm) || cutoff<=min(x, na.rm=na.rm)) {
		stop("cutoff level is not well within range!")
	}
	if(high){
		xl<- x>=cutoff
	} else {
		xl<- x<=cutoff
	}
	
	if(all(xl) | all(!xl)) {
		warnings("nothing to be done, all above/below cutoff !")
		intStart<-intEnd<-integer(0)
		return(list(IS=intStart, IE=intEnd, CUTOFF=cutoff, MG=maxgap, MR=minrun, SPL=splitLen))		
	}
	n<-length(x)
	
	if(!is.null(xPos) && is.null(xRange)){
		if(!is.numeric(xPos)) stop("xPos is not numeric!")
		if(n!=length(xPos)) stop("the length of xPos doesn't match the length of x!")
		if(all((xPos[-1]-xPos[-n])>maxgap) | (xPos[n]-xPos[1])<minrun)
			warnings("maxgap and minrun values are not appropriate !!")
	} else if (is.null(xPos) && !is.null(xRange)){
		if(class(xRange) !='IRanges') stop("xRange is not a IRange object!")
		if(n!=length(xRange)) stop("the length of xRange doesn't match the length of x!")
		if(all(start(xRange)[-1]-end(xRange)[-n] > maxgap) | (end(xRange)[n]- start(xRange)[1]) < minrun)
			warnings("maxgap and minrun values are not appropriate !!")
	}
	
	
	# this version of the maxgap minrun utilize Rle
	# every 2 consecutive runs will have different sign, a merit
	xlrle<-Rle(xl)
	i<-min(which(runValue(xlrle)))
	while(i < max(which(runValue(xlrle)))){

		rv<-runValue(xlrle)
		rl<-runLength(xlrle)	
		erl<-cumsum(rl) #idx of the end of each run
		#srl<-erl-rl+1 # idx of the start of each run, not quite used...
		nr<-length(runValue(xlrle))
		if(i == nr){
			break
		} else if(	(is.null(xPos) && is.null(xRange) && rl[i+1]<=maxgap) ||
						(!is.null(xPos) && is.null(xRange) && xPos[erl[i+1]]-xPos[erl[i]]<=maxgap) ||
						(is.null(xPos) && !is.null(xRange) && start(xRange)[erl[i+1]]-end(xRange)[erl[i]]<=maxgap) && (i+1)<nr ){
			# the next F interval is shorter than maxgap, and there is one T interval still hanging 
			# the gap distance is calculated as the distance between the end of the last True in run i to the last FALSE in the run i+1, half open
			if(	(is.null(xPos) && is.null(xRange) && sum(rl[i:min(i+2, nr)])<=splitLen) ||
				(!is.null(xPos) && is.null(xRange) && xPos[erl[min(i+2, nr)]]-xPos[erl[i]] <=splitLen) ||
				(is.null(xPos) && !is.null(xRange) && end(xRange)[erl[min(i+2, nr)]]-start(xRange)[erl[i]]<=splitLen) ){ 
				#the accumulative length of this interval including this F interval is shorter than the splitLen, merge
				runValue(xlrle)[(i+1)]<-TRUE
			} else {
				#move on and ignore this F interval, no merging
				i<-i+2
			}
		} else {
			i<-i+2
		}
	}
	
	rv<-runValue(xlrle)
	rl<-runLength(xlrle)
	erl<-cumsum(rl) 
	srl<-erl-rl+1 
	nr<-length(runValue(xlrle))
	if(!is.null(xPos)){
		ri<-which(rv & xPos[erl]-xPos[srl]>= minrun)
	} else if (!is.null(xRange)){
		ri<-which(rv & end(xRange)[erl]-start(xRange)[srl]>= minrun)
	} else {
		ri<-which(rv & rl>=minrun)
	}
	if(length(ri)>0){
		intStart<-sapply(ri, function(z) sum(rl[seq_len(z-1)])+1)
		intEnd<-sapply(ri, function(z) sum(rl[seq_len(z)]))
		
		## solve neighbouring probes which are both true, but too far away
		tmp<-splitFarNeighbour(intStart=intStart, intEnd=intEnd, xPos=xPos, xRange=xRange, maxgap=maxgap, minrun=minrun)
		intStart<-tmp$IS
		intEnd<-tmp$IE
	} else {
		intStart<-intEnd<-integer(0)
	}
	return(list(IS=intStart, IE=intEnd, CUTOFF=cutoff, MG=maxgap, MR=minrun, SPL=splitLen))		
}

##################################################
# helper function to split far away neighbouring items
##################################################
splitFarNeighbour<-function(intStart=NULL, intEnd=NULL, xPos=NULL, xRange=NULL, maxgap=Inf, minrun=1){
	if(!is.null(xRange) && class(xRange)=='GRanges' || class(xRange)== 'IRanges') xPos<-NULL
	if(!is.null(xPos) || !is.null(xRange) && !is.null(intStart) && !is.null(intEnd) && length(intEnd)==length(intStart) && !is.infinite(maxgap)){
		if(!is.null(xPos)){
			n<-length(xPos)
			gapStart<-which((xPos[-1]-xPos[-n])>maxgap)
		} else {
			n<-length(xRange)
			gapStart<-which((start(xRange)[-1]-end(xRange)[-n])>maxgap)	
		}
		gapEnd<-gapStart+1
		gapir<-IRanges(start=gapStart, end=gapEnd)
		gapir<-reduce(gapir)
		intN<-length(intStart)
		todel<-numeric(0)
		toadd<-numeric(0)
		for(i in 1:intN){
			splitgap<-which(intStart[i]<=gapStart & intEnd[i]>=gapEnd)
			ngap<-length(splitgap)
			if(ngap>0){
				todel<-c(todel, i)
				for(g in 1:ngap){
					if(ngap>1){
						#there are more than one
						if(g==1){
							# the 1st gap in this interval
							if(	(!is.null(xPos) && xPos[gapStart[splitgap[g]]]-xPos[intStart[i]]>=minrun) ||
								(!is.null(xRange) && end(xRange)[gapStart[splitgap[g]]]-start(xRange)[intStart[i]]>=minrun) ){
								#split first half	
								toadd<-c(toadd, intStart[i], gapStart[splitgap[g]])
							}
						}
						if (g==ngap){
							# the last gap in this interval
							if(	(!is.null(xPos) && xPos[intEnd[i]]-xPos[gapEnd[splitgap[g]]]>=minrun) ||
								(!is.null(xRange) && end(xRange)[intEnd[i]]-start(xRange)[gapEnd[splitgap[g]]]>=minrun) ){
								#split second half
								toadd<-c(toadd, gapEnd[splitgap[g]], intEnd[i])
							}
						} else{
							if(	(g+1)<ngap &&	(!is.null(xPos) && xPos[gapStart[splitgap[g+1]]]-xPos[gapEnd[splitgap[g]]]>=minrun) ||
								(!is.null(xRange) && end(xRange)[gapStart[splitgap[g+1]]]-start(xRange)[gapEnd[splitgap[g]]]>=minrun) ){
								#split second half
								toadd<-c(toadd, gapEnd[splitgap[g]], gapStart[splitgap[g+1]])
							}
						}
					} else {
						# there is only one gap within this interval
						if(	(!is.null(xPos) && xPos[gapStart[g]]-xPos[intStart[i]]>=minrun) ||
							(!is.null(xRange) && end(xRange)[gapStart[g]]-start(xRange)[intStart[i]]>=minrun) ){
							#split first half	
							toadd<-c(toadd, intStart[i], gapStart[splitgap[g]])
						}
						if( (!is.null(xPos) && xPos[intEnd[i]]-xPos[gapEnd[splitgap[g]]]>=minrun) ||
							(!is.null(xRange) && end(xRange)[intEnd[i]]-start(xRange)[gapEnd[splitgap[g]]]>=minrun) ){
							#split second half
							toadd<-c(toadd, gapEnd[splitgap[g]], intEnd[i])
						}
					}							
				}
			}
		}
		if(length(todel)>0){
			intStart<-intStart[-todel]
			intEnd<-intEnd[-todel]
		}
		if(length(toadd)>0){
			intStart<-c(intStart, toadd[seq(1, length(toadd), 2)])
			intEnd<-c(intEnd, toadd[seq(2, length(toadd), 2)])
		}	
	}
	return(list(IS=intStart, IE=intEnd))
}

gammaFit<-function(x, wt=NULL){
	if(is.null(wt)) wt <- rep(1,length(x))
	if(length(x) != length(wt)) stop("length of x and wt differ!")
	
	tmp <- cov.wt(data.frame(x),wt=wt)
    shape <- as.numeric(tmp$center/sqrt(tmp$cov))^2
    scale<-as.numeric(tmp$center)/shape
#	par<-optim(c(shape,scale),function(par) sum(dgamma(x,shape=par[1],scale=par[2], log=TRUE)*wt, na.rm=TRUE)  ,lower=.Machine$double.eps, method='L-BFGS-B', control=list(fnscale=-1,maxit=1000))$par
#	return(c(shape=par[1], scale=par[2]))
	return(c(shape=shape, scale=scale))
}

poisFit <- function(x, wt=NULL, maxshift=0) {  	
	if(maxshift>min(x)) stop("maxshift can't be greater than the minimum of x !")
	if(is.null(wt)) wt <- rep(1/length(x),length(x))
	if(length(x) != length(wt)) stop("length of x and wt differ!")
	
	shift<-which.max(sapply(0:maxshift, function(s) sum(dpois(x = x-s, lambda=sum((x-s) * wt, na.rm=T) ,log=TRUE) * wt, na.rm=T)))-1
	lambda <- sum((x-shift) * wt, na.rm=T) 
    return(c(lambda=lambda, shift=shift))
}

nbinomFit <- function(x, wt=NULL, maxshift=0) {  	
	if(maxshift>min(x)) stop("maxshift can't be greater than the minimum of x !")
	if(is.null(wt)) wt <- rep(1,length(x))
	if(length(x) != length(wt)) stop("length of x and wt differ!")
	
	shift <- which.max(sapply(0:maxshift,  function(s) nbinomCLLDD(x, wt, s)$value))-1
    res<-c(nbinomCLLDD(x, wt, shift)$par, shift)
    names(res)<-c('size', 'mu', 'shift')
    return(res)
}


nbinomCLLDD<-function(x, wt=NULL, s=0){
	if(is.null(wt)) wt <- rep(1,length(x))
	if(length(x) != length(wt)) stop("length of x and wt differ!")
	tmp <- cov.wt(matrix(as.vector(x-s)),wt=wt)
	m <- as.numeric(tmp$center)
	v <- as.numeric(tmp$cov)
	size <- if (v > m) m^2/(v - m) else 100
	value<-sum(dnbinom(x = as.vector(x-s), size=size, mu=m ,log=TRUE) * wt, na.rm=T)
	return(list(value=value, par=c(size, m)))
#	optim(c(size,m),function(par) sum(dnbinom(x-s,size=par[1],mu=par[2],log=TRUE)*wt, na.rm=TRUE) , lower=.Machine$double.eps, method='L-BFGS-B', control=list(fnscale=-1,maxit=1000))
}

tmvtfFit<-function(x, wt=NULL){
	#Aeschliman, C., Park, J., & Kak, A. (2010). A novel parameter estimation algorithm for the multivariate t-distribution and its application to computer vision. Computer Vision-ECCV 2010, 594-607.
	if(!is.matrix(x))	x<-matrix(x)
	n<-nrow(x)
	p<-ncol(x)
	if(is.null(wt)) wt <- rep(1,n)
	if(n != length(wt)) stop("rows of x and length of wt differ!")
	tmp <- cov.wt(x,wt=wt)
	c<-unname(tmp$center)
	z<-apply(x, 1, function(r) log(sum((r - c)^2)))
	zbar<-mean(z)
	b<-sum((z-zbar)^2)/n - trigamma(p/2)
	vhat<-floor((1+sqrt(1+4*b))/b)
	return(list(mu=c, df=ifelse(vhat>1, vhat, 1), var=unname(tmp$cov)))
}

avgFunc<-function(x, avg.m='median', trim=0, na.rm=TRUE){
	switch(avg.m,
		mean=mean(x, trim=trim, na.rm=na.rm),
		median=median(x, na.rm=na.rm),
		stop("invalid 'avg.m' specified!")
	)
}
