biomvRGviz<-function(exprgr, gmgr=NULL, prange=NULL, regionID=NULL, seggr=NULL, plotstrand='+', eps=TRUE, tofile=TRUE, ...){
	# exprgr, probe info, with first data column as expression value
	# gmgr,  related annotation data, optional, a TYPE mcol has to be there...fragile
	# seggr, segmentation info, optional, a STATE mcol has to be there...fragile
	# prange, a range to plot, optional, a reasonable region is strongly advisable.
	#Ideogram track obviously is not avaliable for pig on UCSC
	#ideoTrack <- IdeogramTrack(genome = "susScr3", chromosome = "chrX")
	options(ucscChromosomeNames=FALSE)
	withstrand<-FALSE
	#check if there are more chrs, and no prange
	if(length(unique(as.character(seqnames(exprgr))))>1 && is.null(prange)){
		stop("More than 1 chr in exprgr, yet no plot region defined! ")
	}
	if(is.null(prange)){
		prange<-c(unique(as.character(seqnames(exprgr))), floor(min(start(exprgr))/1000)*1000, ceiling(max(end(exprgr))/1000)*1000)
	}
	if(! plotstrand %in% c('+', '-', '*')) stop("Invalid plotstrand parameter specified, must be one of '+' / '-' / '*' !")
	# handle the colour automatically according to gmgr and seggr
	colors<-c('cyan', 'tomato', 'green','purple','gold', 'violet')
	typecode<-NULL
	if(!is.null(gmgr)){
		if(length(unique(values(gmgr)[,'TYPE'])) > length(colors)) stop("There are too many unique levels in values(gmgr)[,'TYPE'], please re-check!")
		typecode<-c(typecode,unique(values(gmgr)[,'TYPE']))
	}
	if(!is.null(seggr)){
		if(length(unique(values(seggr)[,'STATE'])) > length(colors)) stop("There are too many unique levels in values(seggr)[,'STATE'], please re-check!")
		typecode<-unique(c(typecode, unique(values(seggr)[,'STATE'])))
	}
	if(!is.null(typecode)){
		typecode<-typecode[order(typecode)]
		params <- as.list(colors[seq_along(typecode)])
		names(params)<-typecode
	} else {
		params<-list()
	}
	
	if(hasArg(ylab))  ylab <- list(...)$ylab else ylab<-NULL
	if(hasArg(main))  main <- list(...)$main else main<-NULL
	if(hasArg(cex))  cex <- list(...)$cex else cex<-1.5
	if(hasArg(width))  width <- list(...)$width else width<-16
	if(hasArg(height))  height <- list(...)$height else height<-9
	if(hasArg(fontsize))  fontsize <- list(...)$fontsize else fontsize<-9
	if(hasArg(showId))  showId <- list(...)$showId else showId<-TRUE
	
	
	regionID<-ifelse(is.null(regionID), '', paste(regionID, '@', sep=''))
	if(is.null(main))	main<-paste(regionID, prange[1], '.', prange[2],'-', prange[3], '@', paste(colnames(mcols(exprgr)), collapse='&'),  sep='')
	trackList<-list()
		
	# datatrack + 
	if(plotstrand == '+' | plotstrand == '*'){
		ylabp<-ifelse(is.null(ylab), paste('   ', sep=''), paste(ylab, ifelse(withstrand, '+', ''), sep=''))
		dpTrack <- DataTrack(exprgr[seqnames(exprgr)==prange[1] & (strand(exprgr)=='+' | strand(exprgr)=='*') & start(exprgr) >= as.numeric(prange[2]) & end(exprgr) <= as.numeric(prange[3])],  
					name = ylabp, background.title = "darkblue", type = c('p'),  legend = TRUE, groups = colnames(mcols(exprgr)))				
		trackList<-append(trackList, dpTrack)
		if(!is.null(seggr)){
			# segmentation + as a separate state annotation
			# no id, with legend
			segp<-seggr[seqnames(seggr)==prange[1] & (strand(seggr)=='+' | strand(seggr)=='*') & start(seggr) >= as.numeric(prange[2]) & end(seggr) <= as.numeric(prange[3])]
			for(sn in unique(mcols(segp)[,'SAMPLE'])){
				sni<-mcols(segp)[,'SAMPLE']==sn
				spTrack<- AnnotationTrack(segp[sni],  name=paste(sn, ifelse(withstrand, '+', ''), sep=''), id=values(segp)[sni,'STATE'] ,background.title = "Gray", background.panel = "#FFFFFF", showFeatureId = showId, shape = "box")
				feature(spTrack)<- values(segp)[sni,'STATE']
				trackList<-append(trackList, spTrack)
			}
		}
		
		if(!is.null(gmgr)){
			# annodat track +
			gmp<-gmgr[seqnames(gmgr)==prange[1] & (strand(gmgr)=='+' | strand(gmgr)=='*') & start(gmgr) >= as.numeric(prange[2]) & end(gmgr) <= as.numeric(prange[3])]
			apTrack<- AnnotationTrack(gmp, group=names(gmp), name=ifelse(length(gmp)==0, '', paste('  ', sep='')), background.title = "brown", background.panel = "#FFFEDB", showId = showId, shape = "box")
			feature(apTrack)<- values(gmp)[,'TYPE'] ## this is a strong requirement
			trackList<-append(trackList, apTrack)
			rm(gmp)
		}				
	}
	
	# axis track
	axisTrack <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE, littleTicks = TRUE)
	trackList<-append(trackList, axisTrack)
	
	if(plotstrand == '-' | plotstrand == '*'){			
		if(!is.null(gmgr)){
			# annodat track -
			gmm<-gmgr[seqnames(gmgr)==prange[1] & (strand(gmgr)=='-' | strand(gmgr)=='*') & start(gmgr) >= as.numeric(prange[2]) & end(gmgr) <= as.numeric(prange[3])]
			amTrack <- AnnotationTrack(gmm, group=names(gmm), name=ifelse(length(gmm)==0, '', paste('  ', sep='')), background.title = "brown", background.panel = "#FFFEDB", showId = showId, shape = "box")
			feature(amTrack)<- values(gmm)[,'TYPE']
			trackList<-append(trackList, amTrack)
			rm(gmm)
		}
		if(!is.null(seggr)){
			# segmentation + as a separate state annotation
			# no id, with legend
			segm<-seggr[seqnames(seggr)==prange[1] &  (strand(seggr)=='-' | strand(seggr)=='*')  & start(seggr) >= as.numeric(prange[2]) & end(seggr) <= as.numeric(prange[3])]
			for(sn in unique(mcols(segm)[,'SAMPLE'])){
				sni<-mcols(segm)[,'SAMPLE']==sn
				smTrack<- AnnotationTrack(segm[sni],  name=paste(sn, ifelse(withstrand, '-', ''), sep=''), id=values(segm)[sni,'STATE'] ,background.title = "Gray", background.panel = "#FFFFFF", showFeatureId = showId, shape = "box")
				feature(smTrack)<- values(segm)[sni,'STATE']
				trackList<-append(trackList, smTrack)
			}
		}
		#datatrack - 
		ylabm<-ifelse(is.null(ylab), paste('   ', sep=''), paste(ylab, ifelse(withstrand, '-', ''), sep=''))
		dmTrack <- DataTrack(exprgr[seqnames(exprgr)==prange[1] & (strand(exprgr)=='-' | strand(exprgr)=='*') & start(exprgr) >= as.numeric(prange[2]) & end(exprgr) <= as.numeric(prange[3])], 
					name = ylabm, background.title = "darkblue", type = c('p'),   legend = TRUE,  groups = colnames(mcols(exprgr)))
		trackList<-append(trackList, dmTrack)	
	}

	# start plotting
	# initial grapic dev.
	
	if(tofile){
		graphics.off()
		if(eps){
			setEPS()
			postscript(paste(main, '.', plotstrand,'.eps', sep=''), paper='special', width=width, height=height, horizontal=F, fonts=c("sans"))
		} else {
			pdf(paste(main, '.', plotstrand, '.pdf', sep=''), width=width, height=height)
		}
		# append general plotting params.
		# min.distance = 0, min.width = 0, to distinguish adjacent feature within the same group, not sure if there will be unexpected side effects or not.
		params <- append(params, list(cex.axis=cex, main=main, cex.main=cex, min.distance = 0, min.width = 0, cex.legend=cex, cex=cex)) 
		do.call(plotTracks,c(list(trackList), params))
		dev.off()
	} else {
		params <- append(params, list(main=main, min.distance = 0, min.width = 0)) 
		do.call(plotTracks,c(list(trackList), params))
		return(c(list(trackList), params)) 
	}
}

