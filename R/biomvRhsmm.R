##################################################
# re-implement HSMM, one more slot to handle distance array
##################################################
biomvRhsmm<-function(x, maxk=NULL, maxbp=NULL, J=3, xPos=NULL, xRange=NULL, usePos='start', emis.type='norm', com.emis=FALSE, xAnno=NULL, soj.type='gamma', q.alpha=0.05, r.var=0.75, useMC=TRUE, cMethod='F-B', maxit=1, maxgap=Inf, tol=1e-06, grp=NULL, cluster.m=NULL, avg.m='median', prior.m = 'cluster', trim=0, na.rm=TRUE){
	## input checking
	# lock.transition / lock.d, lock transition and sojourn #fixme
	# est.method=c('viterbi', 'smooth')
	
	message('Checking input ...')
	if (!is.numeric(x) &&  !is.matrix(x) && !is(x, "GRanges"))
        stop("'x' must be a numeric vector or matrix or a GRanges object.")
    if(is(x, "GRanges")) {
    	xid<-names(values(x))
    	xRange<-x
    	mcols(xRange)<-NULL
    	x<-as.matrix(values(x))
    	if( any(sapply(as.character(unique(seqnames(xRange))), function(s) length(unique(strand(xRange)[seqnames(xRange)==s])))!=1) )
    		stop('For some sequence, there are data appear on both strands !')
    } else if(length(dim(x))==2){
		xid<-colnames(x)
		x<-as.matrix(x)
	} else {
		message('No dim attributes, coercing x to a matrix with 1 column.')
		x <- matrix(as.numeric(x), ncol=1)
		xid<-NULL
	}
	
	# check the case when Rle in x mcols
	if(!is.null(xRange) & length(xRange) != nrow(x)){
		nr<-length(xRange)
		spRle<-TRUE
		if(nr > .Machine$integer.max) stop("Length of input is longer than .Machine$integer.max")
	} else {
		nr<-nrow(x)
		spRle<-FALSE
	}

	nc<-ncol(x)
	if(is.null(xid)){
		xid<-paste('S', seq_len(nc), sep='')
		colnames(x)<-xid
	}
	
	if (is.null(cMethod) || !(cMethod %in% c('F-B', 'Viterbi'))) 
		stop("'cMethod' must be specified, must be either 'F-B' or 'Viterbi'!")
		
	if (is.null(prior.m) || !(prior.m %in% c('quantile', 'cluster'))) 
		stop("'prior.m' must be specified, must be either 'quantile' or 'cluster'!")	

	if (is.null(emis.type) || !(emis.type %in% c('norm', 'mvnorm', 'pois', 'nbinom', 'mvt', 't'))) 
		stop("'emis.type' must be specified, must be one of 'norm', 'mvnorm', 'pois', 'nbinom', 'mvt', 't'!")
		
	if (prior.m == 'cluster'){
		if(length(find.package('cluster', quiet=T))==0) {
			warning("'cluster' is not found, fall back to quantile method!!!")
			prior.m <- 'quantile'
		} else {
			require(cluster)
		}
	}

	## some checking on xPos and xRange, xRange exist then xPos derived from xRange,
	if(!is.null(xRange) && (is(xRange, "GRanges") || is(xRange, "IRanges")) && !is.null(usePos) && length(xRange)==nr && usePos %in% c('start', 'end', 'mid')){
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
	if (!is.null(maxbp) && (!is.numeric(maxbp) || (length(maxbp) != 1) || (maxbp <= 1) ||  ( !is.null(xPos) && maxbp > max(xPos,na.rm=na.rm)-min(xPos, na.rm=na.rm)))) 
	 	 stop(sprintf("'maxbp' must be a single integer between 2 and the maximum length of the region if xPos is available!"))	
	
	# check grp setting, cluster if needed, otherwise treat as one group
	if(!is.null(grp)) grp<-as.character(grp)
	grp<-preClustGrp(x, grp=grp, cluster.m=cluster.m) #?todo there could be problem if spRle and need clustering
	
	# initial sojourn setup unify parameter input / density input,  using extra distance, non-integer value can give a dtype value
	if(!is.null(xAnno) && !is.null(soj.type) && soj.type %in% c('gamma', 'pois', 'nbinom') && (is(xAnno, "TxDb") || is(xAnno, "GRanges") || is(xAnno, "GRangesList") || is.list(xAnno))){
		#	this is only used when the xAnno object contains appropriate annotation information which could be used as prior for the sojourn dist in the new HSMM model
		# if xAnno is also present, then J will be estimated from xAnno, and pop a warning, ## this only make sense if difference exist in the distribution of sojourn of states.	
		# a further list object allow direct custom input for initial sojourn dist parameters., e.g. list(lambda=c(10, 50, 1000))
		soj<-sojournAnno(xAnno, soj.type=soj.type)
		J<-soj$J
		message('Estimated state number from xAnno: J = ', J)	
		#now, if there is xAnno, then J, maxbp could be inferred, and if xPos exists, then maxk could also be inferred.
		if(is.null(maxbp)){
			# estimating a reasonable number for maxbp
			maxbp<- switch(soj.type,
				nbinom = ceiling(median(soj$mu)),
				pois = ceiling(median(soj$lambda)),
				gamma = ceiling(median(soj$shape*soj$scale)),
				stop("Invalid argument value for 'soj.type'! ")
			)
		}
#		soj<-append(soj, maxbp=maxbp) 
	} else if(is.numeric(J) && J>1){
		message('xAnno is not present or not supported, try to use maxbp/maxk in the uniform prior for the sojourn distribution.') 
		# J is ok
		if(is.null(xPos)){
			# no position as well, in turn means no xRange nor multiple seq, 
			# init if J and maxk are ok
			if (!is.null(maxk) && is.numeric(maxk) && (length(maxk) == 1) && (maxk > 1) &&  (maxk < nr)) {
				message('maxbp and xPos are not present or not valid, using maxk for the sojourn distribution.')
				soj<-list(type = soj.type, J=J, maxk=maxk)
			} else {
				stop(sprintf("'maxk' must be a single integer between 2 and the number of rows of 'x': %d.", nr))
			}
		} else if(!is.null(maxbp) && maxbp > 1){
			# has position and good maxbp, will init it latter, maxk will be estimated there
			soj<-list(J=J, maxbp=maxbp, type = soj.type)
		} else if (!is.null(maxk) && is.numeric(maxk) && (length(maxk) == 1) && (maxk > 1) &&  (maxk < nr)) {
			#has pos, but no good maxbp
				warning('Has positions but no maxbp, using maxk for the sojourn distribution !!!')
				soj<-list(type = soj.type, J=J, maxk=maxk)
		} else {
			stop(sprintf("Both maxk and maxbp are not available!"))
		}
	} else {
		# no good J
		stop("J must be specified or estimated from xAnno !!")
	}
	# so far, soj is a list object, depending on which case
	# case1, soj param J, maxbp from xAnno
	# case2, no pos, has input maxk and J and d, ready for initSojDd # may reinitialize maxk later
	# case3, has input xPos and maxbp

	# check mv vs iterative
	 if (emis.type=='mvnorm' || emis.type=='mvt' ){
		iterative<-FALSE
	} else {
		iterative<-TRUE
	}
	
	## build xRange if not a GRanges for the returning object
	if(is.null(xRange) || !is(xRange, "GRanges")){
		if(!is.null(xRange) && is(xRange, "IRanges")){
			xRange<-GRanges(seqnames='sampleseq', xRange)	
		} else 	if(!is.null(xPos)){
			xRange<-GRanges(seqnames='sampleseq', IRanges(start=xPos, width=1))	
		} else {
			xRange<-GRanges(seqnames='sampleseq', IRanges(start=seq_len(nr), width=1))	
		}
	}
	# get seqnames status	
	seqs<-unique(as.character(seqnames(xRange))) # speed gain for large no. of contig

	### pre calc emis if common prior for all seq, done outside of seqname loop
	message("Estimating common emission prior ...")
	if(com.emis){
		if(iterative){
			if(spRle){
				emis<-lapply(seq_len(ncol(x)), function(c) estEmis(x[,c][[1]], J=J, emis.type=emis.type, prior.m=prior.m, q.alpha=q.alpha, r.var=r.var)) 
			} else {
				emis<-lapply(seq_len(ncol(x)), function(c) estEmis(x[,c], J=J, emis.type=emis.type, prior.m=prior.m, q.alpha=q.alpha, r.var=r.var)) 
			}
		} else {
			if(spRle){
				stop('Rle like structure is not currently supported for mvnorm / mvt!')
			} else {
				emis<-lapply(unique(grp), function(g) estEmis(x[,grp==g], J=J, emis.type=emis.type, prior.m=prior.m, q.alpha=q.alpha, r.var=r.var)) 
			}
			names(emis)<-unique(grp)
		}
	}

	#make it parallel
	if(useMC & length(find.package('parallel', quiet=T))==0) {
		warning("'parallel' is not found, use normal 'lapply' function!!!")
		useMC<-FALSE
		mylapply<-lapply
	} else if (useMC){
		require(parallel)
		mylapply<-mclapply
	} else {
		mylapply<-lapply
	}
	# we have more than one seq to batch
	mcres<-mylapply(seq_along(seqs), function(s) {
		r<-which(as.character(seqnames(xRange)) == seqs[s])
			
		if(length(r)<2){
			warning('Region too short for seq ', s, ' skipped!')
			if(iterative){
				runout<- sapply(seq_len(nc), function(c) list(NA))
			} else {
				runout<- sapply(unique(grp), function(g) list(NA))
			}
		} else {
			runout<-list()
			# prep soj for the c loop, since there are multiple seq, which also means there must be xpos and maxbp
			message(sprintf("Preparing sojourn prior for seq  '%s' ...", seqs[s]))				
			if(is.null(soj$maxk)){
				# either has soj parameter, or has pos and maxbp for unif
				ssoj<-append(soj, initDposV(xPos[r], maxbp))
				if(is.null(ssoj$fttypes)){
					# dont't have soj param
					ssoj<-append(ssoj, list(d=unifMJ(ssoj$maxk*length(r), J)))
				}
			} else {
				# maxk and unif d, no maxbp or xAnno, using the min of input maxk and length of the seq
				ssoj<-append(soj, list(d=unifMJ(min(soj$maxk,length(r)), J)))
				ssoj$maxk<-min(soj$maxk,length(r))
			}	
			ssoj <- initSojDd(ssoj)
	
			for(g in unique(grp)){
				gi<-grp==g
				if(iterative){
					for(ci in which(gi)){
						message(sprintf("Building HSMM for seq '%s' in column '%s' ...", seqs[s], ci))				
						if(com.emis) {
							semis<-emis[[ci]]
						}else {
							if(spRle){
								semis<-estEmis(x[,ci][[1]][r], J=J, prior.m=prior.m, emis.type=emis.type, q.alpha=q.alpha, r.var=r.var)
							} else {
								semis<-estEmis(x[r,ci], J=J, prior.m=prior.m, emis.type=emis.type, q.alpha=q.alpha, r.var=r.var)
							}
						}
						if(spRle){
							grunout<-tryCatch(hsmmRun(x[,ci][[1]][r], xid[ci], xRange[r], ssoj, semis, cMethod, maxit, maxgap,  tol, avg.m=avg.m, trim=trim, na.rm=na.rm, com.emis=com.emis), error=function(e){ return(e) })
						} else {
							grunout<-tryCatch(hsmmRun(x[r,ci], xid[ci], xRange[r], ssoj, semis, cMethod, maxit, maxgap,  tol, avg.m=avg.m, trim=trim, na.rm=na.rm, com.emis=com.emis), error=function(e){ return(e) }) 	
						}
						if(length(grep('error', class(grunout), ignore.case=T))!=0){
							warning('Model failed for seq ', s, ' and column ', c, '.\n', grunout)
							runout<-append(runout, list(NA))
						} else {
							runout<-append(runout, list(grunout)	)
						}	
					}
				} else {
					message(sprintf("Building HSMM for group '%s' ...", g))
					if(com.emis) {
						semis<-emis[[g]]
					}else {
						if(spRle){
							stop('Rle like structure is not currently supported for mvnorm / mvt!')
						} else {
							semis<-estEmis(x[r,gi], J=J, prior.m=prior.m, emis.type=emis.type, q.alpha=q.alpha, r.var=r.var)
						}
					}
					grunout<-tryCatch(hsmmRun(x[r,gi], xid[gi], xRange[r], ssoj, semis, cMethod, maxit, maxgap, tol, avg.m=avg.m, trim=trim, na.rm=na.rm, com.emis=com.emis)	, error=function(e){ return(e) })
					if(length(grep('error', class(grunout), ignore.case=T))!=0){
						warning('Model failed for seq ', s, ' and group ', g, '.\n', grunout)
						runout<-append(runout, list(NA))
					} else {
						runout<-append(runout, list(grunout))
					}
				}
			} # end for g
		} # end if 
		return(runout)
	})
	
	message("Building HSMM complete, preparing output ...")
	res<-do.call(c, lapply(seq_along(seqs), function(s){
		do.call(c, lapply(seq_len(ifelse(iterative, nc, length(unique(grp)))), function(g){
			 if(!is.na(mcres[[s]][[g]][1])){
			 	return(mcres[[s]][[g]]$res)
			 } else {
			 	return(GRanges())
			 }
		}))
	}))
	
	seqr<-table(seqnames(xRange))
	ssp<-do.call(rbind, lapply(seq_along(seqs), function(s){
		do.call(cbind, lapply(seq_len(ifelse(iterative, nc, length(unique(grp)))), function(g){
			 if(!is.na(mcres[[s]][[g]][1])){
			 	return(DataFrame(s=mcres[[s]][[g]]$yhat, sp=mcres[[s]][[g]]$yp))
			 } else {
			 	return(DataFrame(s=Rle(NA, seqr[seqs[s]]), sp=Rle(NA, seqr[seqs[s]])))
			 }
		}))
	}))
	
	emis.par<-do.call(rbind, lapply(seq_along(seqs), function(s){
		do.call(c, lapply(seq_len(ifelse(iterative, nc, length(unique(grp)))), function(g){
			 if(!is.na(mcres[[s]][[g]][1])){
			 	return(list(mcres[[s]][[g]]$emispar))
			 } else {
			 	return(list(NA))
			 }
		}))
	}))
	rownames(emis.par)<-seqs
	colnames(emis.par)<-if(iterative) xid else unique(grp)
	
	soj.par<-do.call(rbind, lapply(seq_along(seqs), function(s){
		do.call(c, lapply(seq_len(ifelse(iterative, nc, length(unique(grp)))), function(g){
			 if(!is.na(mcres[[s]][[g]][1])){
			 	return(list(mcres[[s]][[g]]$sojpar))
			 } else {
			 	return(list(NA))
			 }
		}))
	}))
	rownames(soj.par)<-seqs
	colnames(soj.par)<-if(iterative) xid else unique(grp)

	# setup input data and state to xRange for returning
	values(xRange)<-DataFrame(do.call(DataFrame, lapply(seq_len(nc), function(c) x[,c])), ssp, row.names = NULL)
	colnames(mcols(xRange)) <- c(xid, paste(rep(c('S', 'SP'), ifelse(iterative, nc, length(unique(grp)))), rep(if(iterative) xid else unique(grp), each=2), sep='.'))
	new("biomvRCNS",  
		x = xRange, res = res,
		param=list(J=J, maxk=maxk, maxbp=maxbp, maxgap=maxgap, soj.type=soj.type, emis.type=emis.type, q.alpha=q.alpha, r.var=r.var, iterative=iterative, cMethod=cMethod, maxit=maxit, tol=tol, grp=grp, cluster.m=cluster.m, avg.m=avg.m, prior.m=prior.m, trim=trim, na.rm=na.rm, soj.par=soj.par, emis.par=emis.par)
	)
}



hsmmRun<-function(x, xid='sampleid', xRange, soj, emis, cMethod='F-B', maxit=1, maxgap=Inf, tol= 1e-6, avg.m='median', trim=0, na.rm=TRUE, com.emis=FALSE){
	# now x should be a matrix, when emis.type=mvt or mvnorm, ncol>1
	if(is.null(dim(x))) x<-matrix(as.vector(x))
	colnames(x)<-xid
	nr<-nrow(x)
	J<-soj$J
	maxk<-soj$maxk
	s<-unique(as.character(seqnames(xRange))) 
		
	# create default uniform initial probability
	init<-rep(1/J, J) # start with uniform
	# create default uniform transition probability
	trans <- matrix(1/(J-1), nrow = J, ncol=J)
	diag(trans)<-0
	
	#estimation of most likely state sequence
	#define likelihood
	ll <- rep(NA,maxit)
	# start MM iteration
	for(it in 1:maxit) {
		message(sprintf("[ hsmmRun ] seq '%s' column '%s' iteration: %d", s, paste(xid, collapse='+'), it))
		# reestimationg of emmision   
		emis<-initEmis(emis=emis, x=x)
		B  = .C("backward", a=as.double(trans), pi=as.double(init), b=as.double(emis$p), d=as.double(soj$d), D=as.double(soj$D),
				  maxk=as.integer(maxk), DL=as.integer(nrow(soj$d)), T=as.integer(nr), J=as.integer(J), 
				  eta = double(nrow(soj$d)*J), L=double(nr*J), N=double(nr), ahat=double(J*J), pihat=double(J),
				  F=double(nr*J), G=double(nr*J), L1 = double(nr*J), si=double(nr*J), PACKAGE='biomvRCNS')

		#update initial prob PI, transition >=0 check
#		init<-B$pihat
		init<- abs(B$pihat)/sum(abs(B$pihat))
		trans <- matrix(B$ahat,ncol=J)
		trans[trans<0] <- 0
		
		#check B$L
		if(all(is.nan(B$L))) {
		  stop("Sojourn distribution does not work well, NaN in B$L ")
		}
		#update emission according to the new estimated distribution parameters using B$L
		if(!com.emis){
				emis<-initEmis(emis=emis, x=x, B=B)
		}
		# update sojourn dD, using B$eta
		soj<-initSojDd(soj=soj, B=B)
		
		# log-likelihood for this it, using B$N
		ll[it]<-sum(log(B$N))
		if( it>1 && abs(ll[it]-ll[it-1]) < tol) {
			break()	
		}
	}	 # end for maxit
	
	BL<-matrix(B$L,ncol=J)
	BL<-t(apply(BL, 1, function(x) abs(x)/sum(abs(x))))
	
	# switch cMethod
	if (cMethod=='Viterbi'){
		emis<-initEmis(emis=emis, x=x)
		
		logtrans<-log(trans)
		logtrans[logtrans==-Inf] <- -.Machine$double.xmax
		loginit<-log(init)
		loginit[loginit==-Inf] <- -.Machine$double.xmax
		logd = log(soj$d)
		logd[logd==-Inf] <- -.Machine$double.xmax
		logD = log(soj$D)
		logD[logD==-Inf] <- -.Machine$double.xmax
		logb<-log(emis$p)
		logb[logb==-Inf] <- -.Machine$double.xmax
		
		V  = .C("logviterbi", a=as.double(logtrans), pi=as.double(loginit), b=as.double(logb), d=as.double(logd), D=as.double(logD),
          maxk=as.integer(maxk), DL=as.integer(nrow(soj$d)), T=as.integer(nr), J=as.integer(J), 
          alpha = double(nr*J), shat=integer(nr), si=double(nr*J), opt=integer(nr*J), ops=integer(nr*J), PACKAGE='biomvRCNS')
        yhat<-V$shat+1      
	} else if (cMethod=='F-B'){
		## assign states and split if necessary.
		yhat<-apply(BL, 1, which.max)
	}
	yp<-Rle(as.vector(BL)[(yhat-1)*nr+1:nr])
	if(!is.null(soj$fttypes)){
		yhat<-soj$fttypes[yhat]
	}
	yhat<-Rle(as.character(yhat))
	
	# setup this new res gr
	Ilist<-lapply(unique(yhat), function(j) do.call(cbind, splitFarNeighbouryhat(yhat, xRange=xRange, maxgap=maxgap, state=j)))
	names(Ilist)<-unique(yhat)
	res<- do.call('c', 
					lapply(unique(yhat), function(j)
						GRanges(seqnames=as.character(seqnames(xRange)[1]), 
							IRanges(start=rep(start(xRange)[Ilist[[j]][,'IS']], length(xid)), end=rep(end(xRange)[Ilist[[j]][,'IE']], length(xid))), 
							strand=rep(strand(xRange)[Ilist[[j]][,'IS']], length(xid)), 
							SAMPLE=Rle(xid, rep(nrow(Ilist[[j]]), length(xid))), 
							STATE=Rle(as.character(j), nrow(Ilist[[j]])*length(xid)), 
							AVG=Rle(as.numeric(sapply(xid, function(s) apply(Ilist[[j]], 1, function(r) apply(as.matrix(x[r[1]:r[2],s]), 2, avgFunc, avg.m=avg.m, trim=trim, na.rm=na.rm)))))
						)
					)
			)
	seqlevels(res)<-seqlevels(xRange)
	
	# return the soj. and emis par, maybe also the estimated state emis$p
	emispar2ret<-c('mu', 'var', 'size', 'df')
	ei<-names(emis) %in% emispar2ret
	sojpar2ret<-c('lambda', 'shift', 'mu', 'size', 'scale', 'shape')
	si<-names(soj) %in% sojpar2ret
			
	return(list(yhat=yhat, res=res, yp=yp, emispar=emis[ei], sojpar=soj[si]))
}


##################################################
#			HSMM helper functions
##################################################
##################################################
# create lists of uniform prior
##################################################
unifMJ<-function(M,J, ints=NULL){
	# not finished create unif with supplied ints
	if(is.null(ints)){
		ret<-do.call(cbind, lapply(1:J, function(j) dunif(1:M, 0, M)))		
	} else if(min(ints)<1 | max(ints)> M) {
		stop('All supplied intervals should be in the range of [1, M]')
	} else {
		# ints, a dataframe marks the starts and ends for each interval
		ret<-do.call(cbind, lapply(1:J, function(j) do.call(c, lapply(1:nrow(ints), function(i) dunif(1:M, ints[i,1], ints[i,2])))))		
	}
	ret
}




##################################################
# estimate sojourn distribution from annotation using gamma
##################################################
sojournAnno<-function(xAnno, soj.type= 'gamma', pbdist=NULL){ 
	# xAnno has to be a Grange / rangedata obj
	# check if xAnno class, txdb or dataframe or rangedata ...
	# there is also the possibility of proposing an empirical number for the states.
	# must ensure there are at least 2 for each state ? todo
		
	if(is(xAnno, "TxDb")) {   
	   if(length(find.package('GenomicFeatures', quiet=T))==0) {
			stop("'GenomicFeatures' is not found !!!")
		} else {
			require(GenomicFeatures)
		}
		J<-3
		fttypes<-c('intergenic', 'intron', 'exon')
		#3 feature type, exon, intron, intergenic
		transc <- transcripts(xAnno) # this give you all cds ranges ungrouped
		intergenic<-gaps(transc)
		
		# gaps() will by default produce extra * ranges and full range for empty chr
		# https://stat.ethz.ch/pipermail/bioconductor/2013-May/052976.html
		intergenic<-intergenic[strand(intergenic)!='*']
		intergenic<-intergenic[which(width(intergenic) != seqlengths(intergenic)[as.character(seqnames(intergenic))])]
		
		exon <- exons(xAnno) # this give you all exon ranges ungroupped
		intron<- unlist(intronsByTranscript(xAnno))
		ftdist<-list(intergenic=width(intergenic), intron=width(intron), exon=width(exon))
	} else if(is(xAnno, "GRanges")){
		# then the first column of the elementMetadata must be a character vector marking the type of features.
		# need a way to sort 
		fts<-values(xAnno)[,1]
		fttypes<- unique(fts)
		J<-length(fttypes)
		ftdist<-lapply(1:J, function(j) width(xAnno[fts==fttypes[j]]))
	} else if (is(xAnno, "GRangesList")){
		ng<-length(xAnno)
		## the assumptions are number of ft types could be different from different list entry, group wise analysis, which of coz could be wrapped using foreach(ng) single group approach
		fts<-lapply(1:ng, function(g) values(xAnno[[g]])[,1])
		fttypesL<- lapply(fts, function(x) unique(x[order(x)]))
		J<-length(unique(unlist(fts)))
		ftdist<-lapply(1:ng, function(g) lapply(1:length(fttypesL[[g]]), function(j) width(xAnno[[g]][fts[[g]]==fttypesL[[g]][j]])))
		fttypes<-unique(unlist(fttypesL))
		fttypes<-fttypes[order(fttypes)]
		ftdist<-lapply(fttypes, function(t) unlist(lapply(1:ng, function(g) ftdist[[g]][[match(t, fttypesL[[g]])]])))
	} else if (is.list(xAnno)){
		# checking if input list has valid specs.
		paramID<- names(xAnno) 
		#check name fall into the right pool
		paramIDok<- switch(soj.type,
			nbinom = all( paramID %in% c('mu', 'size', 'shift')),
			pois = all( paramID %in% c('lambda', 'shift')),
			gamma = all( paramID %in% c('scale', 'shape')),
			stop("Invalid argument value for 'soj.type'! ")
		)
		#check length of vectors match
		J<-unique(sapply(xAnno, length))
		if(length(J)!=1) {
			stop("Length of vectors in xAnno are not equal, can't get valid state number J!")
		}
		soj<-append(xAnno, list(type = soj.type, J=J))
		return(soj)
	}
	
	# should switch here between different soj.type	
	soj<-list(type = soj.type, fttypes=fttypes, J=J)
	if(soj.type=='gamma'){
		shape<-numeric()
		scale<-numeric()
		for(j in 1:J){
			param<-gammaFit(ftdist[[j]])
			if(! is.null(pbdist)){
				# if distance between points are even
				param['scale'] <- param['scale'] / pbdist 
			}
			shape<-c(shape, param['shape'])
			scale<-c(scale, param['scale'])
		}
		soj<-append(soj, list(shape=unname(shape), scale=unname(scale)) )
	} else if (soj.type == 'nbinom'){
		size<-numeric()
		mu<-numeric()
		shift<-numeric()
		for(j in 1:J){
			param<-nbinomFit(ftdist[[j]])
			if(! is.null(pbdist)){
				# if distance between points are even
				param['mu'] <-  param['mu'] / pbdist
			}
			size<-c(size, param['size'])
			mu<-c(mu, param['mu'])
			shift<-c(shift, param['shift'])
		}
		soj<-append(soj, list(size=unname(size), mu=unname(mu), shift=unname(shift)) )			
	} else if (soj.type == 'pois'){
		lambda<-numeric()
		shift<-numeric()
		for(j in 1:J){
			param<-poisFit(ftdist[[j]])
			if(! is.null(pbdist)){
				# if distance between points are even
				param['lambda'] <- param['lambda'] / pbdist 
			}
			lambda<-c(lambda, param['lambda'])
			shift<-c(shift, param['shift'])
		}
		soj<-append(soj, list(lambda=unname(lambda), shift=unname(shift)) )			
	}
	#return soj object
	return(soj)
}

initDposV<-function(xpos, maxbp){
	# for each position, find the maxk
	nr<-length(xpos)
	maxbpidx<-sapply(1:nr, function(i) max(which(xpos[i]+maxbp > xpos))) # >= will cause a NA at the end of each position,
	# find the maxk idx
	maxk<-max(maxbpidx - seq_len(nr))+1
	# initialize the maxk position list for each position
	## option 2, a TM * J matrix
	dposV<-c(sapply(1:nr, function(t) xpos[t:(t+maxk-1)]-xpos[t]))+1 # dposV[(t-1)*maxk+u]
	# sub > maxbp and NA
	dposV[which(dposV>maxbp)]<-NA
	return(list(dposV=dposV, maxk=maxk))
}



initSojDd <- function(soj, B=NULL) {
	# take the initial soj dist parameter / or sample of density
	if(! soj$type %in% c('nparam', 'gamma', 'pois', 'nbinom')) stop("invalid sojourn type found in soj$type !")
	# parameter initialisation
	J<-soj$J
	maxk<-soj$maxk
	if(is.null(soj$dposV)){
		dposV<-1:maxk
	} else {
		## here d for all positions should be aggregated within each J
		dposV<-soj$dposV
	}
	idx<- !is.na(dposV)
	dposV[!idx]<- .Machine$integer.max
	nb <- length(dposV)/maxk
	if(soj$type == "gamma") {
		if(!is.null(B)){
			# then this is an update run
			soj$d <- matrix(B$eta+.Machine$double.eps,ncol=J)
			soj$shape <- soj$scale <- numeric(J)
			for(j in 1:J) {           
				param <- gammaFit(dposV[idx],wt=soj$d[idx,j])
				soj$shape[j] <- param['shape']
				soj$scale[j] <- param['scale'] 	  
			}
		}	          
		if(!is.null(soj$shape) && !is.null(soj$scale) && length(soj$shape)==J && length(soj$scale)==J) {
			# for update with para estimated from B, or initial sojourn using param
			soj$d<-sapply(1:J, function(j) dgamma(dposV, shape=soj$shape[j], scale=soj$scale[j]))
		} # else assume soj$d exist.
	} else if (soj$type == "pois") {
		if(!is.null(B)){
			# then this is an update run, re-estimation of dist params
			soj$d <- matrix(B$eta+.Machine$double.eps,ncol=J)
			soj$shift <- soj$lambda <- numeric(J)
			
			ftidx<-apply(soj$d, 2, function(xc) xc >.Machine$double.eps & idx)
			maxshift<- sapply( 
				lapply(1:J, function(j)  which(ftidx[,j]) %% maxk ), 
				function(ic) min(dposV[seq(from=min(ic[ic>0]), to=nb*maxk, by=maxk)])
			)
#			cat(maxshift, '\n')
			for(j in 1:J) { 
				param <- poisFit(dposV[ftidx[,j]], wt=soj$d[ftidx[,j],j], maxshift=maxshift[j])
				soj$lambda[j] <- param['lambda']
				soj$shift[j] <- param['shift']
			}
		}
		if(!is.null(soj$shift) && !is.null(soj$lambda) && length(soj$lambda)==J && length(soj$shift)==J) {
			# for update with para estimated from B, or initial sojourn using param
			soj$d<-sapply(1:J, function(j) dpois(dposV-soj$shift[j], lambda=soj$lambda[j]))
		} # else assume soj$d exist.
		
	}  else if (soj$type == "nbinom") {
		if(!is.null(B)){
			# then this is a update run, re-estimation of dist params
			soj$d <- matrix(B$eta+.Machine$double.eps,ncol=J)
			soj$shift <- soj$size <- soj$mu <- numeric(J)    
			
			ftidx<-apply(soj$d, 2, function(xc) xc >.Machine$double.eps & idx)
			maxshift<- sapply( 
				lapply(1:J, function(j)  which(ftidx[,j]) %% maxk ), 
				function(ic) min(dposV[seq(from=min(ic[ic>0]), to=nb*maxk, by=maxk)])
			)
			for(j in 1:J) { 
				param <- nbinomFit(dposV[ftidx[,j]],wt=soj$d[ftidx[,j],j], maxshift=maxshift[j])
				soj$size[j] <- param['size']
				soj$mu[j] <- param['mu']
				soj$shift[j] <- param['shift']
			}	
		}
		if(!is.null(soj$shift) && !is.null(soj$size) && !is.null(soj$mu) && length(soj$mu)==J && length(soj$shift)==J && length(soj$size)==J) {
			# for update with para estimated from B, or initial sojourn using param
			soj$d<-sapply(1:J, function(j) dnbinom(dposV-soj$shift[j],size=soj$size[j],mu=soj$mu[j]) )
		} # else assume soj$d exist.
	} else if (soj$type == "nparam") {
		if(!is.null(B)){
			 soj$d <- matrix(B$eta+.Machine$double.eps, ncol=J)
		} # else assume soj$d exist.
	}
	soj$d<-sapply(1:J, function(j) sapply(1:nb, function(t) soj$d[((t-1)*maxk+1):(t*maxk),j]/sum(soj$d[((t-1)*maxk+1):(t*maxk),j], na.rm=T)))
	soj$d[is.na(soj$d)]<-0 # NaN enters when all 0 in the den when normalizing 
	# add D slot
	soj$D <- sapply(1:J, function(j) sapply(1:nb, function(t) rev(cumsum(rev(soj$d[((t-1)*maxk+1):(t*maxk),j])))))
	return(soj)
}


##################################################
# to estimate emission par
##################################################
estEmis<-function(x, J=3, prior.m='quantile', emis.type='norm', q.alpha=0.05, r.var=0.75){
	if(is.null(dim(x))) x<-matrix(as.vector(x))
	emis<-list(type=emis.type)
	if (prior.m == 'cluster'){
		xclust<-clara(x, J)
		if(ncol(x)>1){
			if(emis$type == 'norm' || emis$type== 'mvnorm' || emis$type == 'mvt') {
				emis$var<-lapply(order(xclust$medoids[,1]), function(j) cov(x[xclust$clustering==j,]))
			} 
			emis$mu<-lapply(order(xclust$medoids[,1]), function(j) xclust$medoids[j,])
		} else {
			if(emis$type == 'norm' || emis$type== 'mvnorm' || emis$type == 'mvt') {
				emis$var<-sapply(order(xclust$medoids), function(j) var(x[xclust$clustering==j,]))
				emis$var[is.na(emis$var)]<-mean(emis$var[!is.na(emis$var)])
			}
			emis$mu<-as.numeric(xclust$medoids)[order(xclust$medoids)]
		}	
	}
	
	# old quantile method
	if(prior.m == 'quantile'){
		emis$mu <- estEmisMu(x, J, q.alpha=q.alpha)
		if(emis$type == 'norm' || emis$type== 'mvnorm' || emis$type == 'mvt') {
			emis$var <- estEmisVar(x, J, r.var=r.var)
		}
	}	
	if (emis$type == 'nbinom'){
		emis$size <- rep(nbinomCLLDD(x)$par[1], J) # common prior
	}
	if (emis$type == 't' || emis$type == 'mvt'){
		emis$df <- rep(1, J) # common prior
	}
	return(emis)
}

##################################################
# to estimate segment wise mean vector/list
##################################################
estEmisMu<- function(x, J, q.alpha=0.05, na.rm=TRUE){
	nc<-ncol(x)
	if(J >1){
		if(is.null(nc) || nc == 1){
			# univariate
			ret<-as.numeric(quantile(x, seq(from=q.alpha, to=(1-q.alpha), length.out=J), na.rm=na.rm))
		} else {
			# multiple
			muv<-t(sapply(1:nc, function(i) as.numeric(quantile(x[,i], 	seq(from=q.alpha, to=(1-q.alpha), length.out=J), na.rm=na.rm))))
			ret<-lapply(apply(muv, 2, list),unlist)
		}
	} else {
		ret<-mean(x, na.rm=na.rm)
	}
	return(ret)
}
##################################################
# to estimate segment wise variance vector / covariance matrix list
##################################################
estEmisVar<-function(x, J=3, na.rm=TRUE, r.var=0.75){
	# r.var is the expected ratio of variance for state 1 and J versus any intermediate states
	# a value larger than 1 tend to give more extreme states;  a value smaller than 1 will decrease the probability of having extreme state, pushing it to the center.
	nc<-ncol(x)
	f.var<-rep(ifelse(r.var>=1, 1, r.var), J)
	if(J%%2 == 1) {
		f.var[(J+1)/2]<-ifelse(r.var>=1, 1/r.var, 1)
	} else if( J>2 ){
		f.var[(J/2):(J/2+1)]<-ifelse(r.var>=1, 1/r.var, 1)
	}	
	if(J >1){
		if(is.null(nc) || nc == 1){
			# univariate
			ret<-var(x, na.rm=na.rm)*f.var
		} else {
			# multiple
			if(na.rm) na.rm<-'complete.obs'
			ret<-lapply(1:J, function(j) cov(x, use=na.rm)*f.var[j])
		}
	} else {
		ret<-var(as.numeric(x), na.rm=na.rm)
	}
	
	return(ret)
}
##################################################
# initialize and update emission probability
##################################################
initEmis<-function(emis, x, B=NULL){
	if(is.null(B)){
		# then this is for p initialization
		J<-length(emis$mu)
		if(emis$type == 'mvnorm') {
			emis$p <- sapply(1:J, function(j) dmvnorm(x, mean = emis$mu[[j]],  sigma = emis$var[[j]])) # here sigma requires cov mat
		} else if (emis$type == 'pois'){
			emis$p <-sapply(1:J, function(j) dpois(x, lambda=emis$mu[j]))
		} else if (emis$type == 'norm'){
			emis$p <- sapply(1:J, function(j) dnorm(x,mean=emis$mu[j], sd=sqrt(emis$var[j]))) # here is sd
		} else if(emis$type == 'nbinom'){
			emis$p <- sapply(1:J, function(j) dnbinom(x,size=emis$size[j], mu=emis$mu[j]))
		} else if(emis$type == 'mvt') {
			emis$p <- sapply(1:J, function(j) dmvt(x, df=emis$df[[j]], delta = emis$mu[[j]],  sigma = cov2cor(emis$var[[j]]), log=FALSE)) # here sigma requires cor mat
		} else if(emis$type == 't') {
			emis$p <- sapply(1:J, function(j) dt(x, df=emis$df[j], ncp = emis$mu[j])) 
		}
		emis$p<-emis$p/rowSums(emis$p) # normalized
	} else {
		# then this is for the re-estimation of emis param
		J<-B$J
		BL<-matrix(B$L,ncol=J)
		BL<-t(apply(BL, 1, function(x) abs(x)/sum(abs(x))))
		if(emis$type == 'mvnorm') {
			isa <-  !apply(is.na(x),1,any) # Find rows with NA's (cov.wt does not like them)
			tmp <- apply(BL, 2, function(cv) cov.wt(x[isa, ], cv[isa])[c('cov', 'center')]) # x is already a matrix 
			emis$mu <- lapply(tmp, function(l) unname(l[['center']]))
			emis$var <- lapply(tmp, function(l) unname(l[['cov']]))
		} else if (emis$type == 'pois'){
			isa <- !is.na(x)
			emis$mu <- apply(BL, 2, function(cv) weighted.mean(matrix(x[isa]), cv[isa]))
		} else if (emis$type == 'norm'){
			isa <- !is.na(x)
			tmp <- apply(BL, 2, function(cv) unlist(cov.wt(matrix(x[isa]), cv[isa])[c('cov', 'center')]))
			emis$mu <- tmp['center',]
			emis$var <- tmp['cov', ]		
		} else if(emis$type == 'nbinom'){
			isa <- !is.na(x)
			tmp <- apply(BL, 2, function(cv) nbinomFit((x[isa]), cv[isa]))
   			emis$mu <- tmp['mu',]
			emis$size <- tmp['size', ]
		} else 	if(emis$type == 'mvt') {
			isa <-  !apply(is.na(x),1,any) 
			tmp <- apply(BL, 2, function(cv) tmvtfFit((x[isa, ]), cv[isa]))
			emis$mu <- lapply(tmp, function(l) l[['mu']])
			emis$df <- lapply(tmp, function(l) l[['df']])
			emis$var <- lapply(tmp, function(l) l[['var']])
		} else if(emis$type == 't') {
			isa <- !is.na(x)
			tmp <- apply(BL, 2, function(cv) tmvtfFit((x[isa]), cv[isa]))
			emis$mu <- sapply(tmp, function(l) l[['mu']])
			emis$df <- sapply(tmp, function(l) l[['df']])
		} 
	}
	return(emis)
}	

splitFarNeighbouryhat<-function(yhatrle, xPos=NULL, xRange=NULL, maxgap=Inf, state=NULL){
#	yhatrle<-Rle(yhat)
	rv<-runValue(yhatrle)
	rl<-runLength(yhatrle)
	ri<-which(rv==as.character(state))
	intStart<-sapply(ri, function(z) sum(rl[seq_len(z-1)])+1)
	intEnd<-sapply(ri, function(z) sum(rl[seq_len(z)]))

	if(length(intStart)>0 && !is.null(maxgap) && !is.null(xPos) || !is.null(xRange)){
		tmp<-splitFarNeighbour(intStart=intStart, intEnd=intEnd, xPos=xPos, xRange=xRange, maxgap=maxgap)
		intStart<-tmp$IS
		intEnd<-tmp$IE
	}
	return(list(IS=intStart, IE=intEnd))
}

