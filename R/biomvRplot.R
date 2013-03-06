biomvRGviz<-function(exprgr, gmgr=NULL, prange=NULL, regionID='regionID', seggr=NULL, plotstrand='+', eps=TRUE, tofile=TRUE,showId=TRUE, ...){
	# exprgr, probe info, with first data column as expression value
	# gmgr,  related annotation data, optional, a TYPE mcol has to be there...fragile
	# seggr, segmentation info, optional, a STATE mcol has to be there...fragile
	# prange, a range to plot, optional, a reasonable region is strongly advisable.
	#Ideogram track obviously is not avaliable for pig on UCSC
	#ideoTrack <- IdeogramTrack(genome = "susScr3", chromosome = "chrX")
	options(ucscChromosomeNames=FALSE)
	#check if there are more chrs, and no prange
	if(length(unique(as.character(seqnames(exprgr))))>1 && is.null(prange)){
		stop("More than 1 chr in the Grange object, yet no plot region defined! ")
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
	
	if(is.null(main))	main<-paste(regionID, '@', prange[1], '.', prange[2],'-', prange[3], '@', colnames(mcols(exprgr)),  sep='')
	trackList<-list()
		
	# datatrack + 
	if(plotstrand == '+' | plotstrand == '*'){
		if(is.null(ylab)) ylab<-paste('E@+', sep='') else ylab<-paste(ylab, '@+', sep='')
		dpTrack <- DataTrack(exprgr[seqnames(exprgr)==prange[1] & (strand(exprgr)=='+' | strand(exprgr)=='*') & start(exprgr) >= as.numeric(prange[2]) & end(exprgr) <= as.numeric(prange[3])],  name = ylab, background.title = "darkblue", type = c('p'),  legend = TRUE)				
		trackList<-append(trackList, dpTrack)
		if(!is.null(seggr)){
			# segmentation + as a separate state annotation
			# no id, with legend
			segp<-seggr[seqnames(seggr)==prange[1] & (strand(seggr)=='+' | strand(seggr)=='*') & start(seggr) >= as.numeric(prange[2]) & end(seggr) <= as.numeric(prange[3])]
			spTrack<- AnnotationTrack(segp,  name=paste('S', '@+', sep=''), id=values(segp)[,'STATE'] ,background.title = "Gray", background.panel = "#FFFFFF", showFeatureId = showId, shape = "box")
			feature(spTrack)<- values(segp)[,'STATE']
			trackList<-append(trackList, spTrack)
		}
		
		if(!is.null(gmgr)){
			# annodat track +
			gmp<-gmgr[seqnames(gmgr)==prange[1] & (strand(gmgr)=='+' | strand(gmgr)=='*') & start(gmgr) >= as.numeric(prange[2]) & end(gmgr) <= as.numeric(prange[3])]
			apTrack<- AnnotationTrack(gmp, group=names(gmp), name=ifelse(length(gmp)==0, '', paste('G', '@+', sep='')), background.title = "brown", background.panel = "#FFFEDB", showId = showId, shape = "box")
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
			amTrack <- AnnotationTrack(gmm, group=names(gmm), name=ifelse(length(gmm)==0, '', paste('G', '@-', sep='')), background.title = "brown", background.panel = "#FFFEDB", showId = showId, shape = "box")
			feature(amTrack)<- values(gmm)[,'TYPE']
			trackList<-append(trackList, amTrack)
			rm(gmm)
		}
		if(!is.null(seggr)){
			# segmentation + as a separate state annotation
			# no id, with legend
			segm<-seggr[seqnames(seggr)==prange[1] &  (strand(seggr)=='-' | strand(seggr)=='*')  & start(seggr) >= as.numeric(prange[2]) & end(seggr) <= as.numeric(prange[3])]
			smTrack<- AnnotationTrack(segm,  name=paste('S', '@-', sep=''),  id=values(segm)[,'STATE'], background.title = "Gray", background.panel = "#FFFFFF", showFeatureId = showId, shape = "box")
			feature(smTrack)<- values(segm)[,'STATE']
			trackList<-append(trackList, smTrack)
		}
		#datatrack - 
		if(is.null(ylab)) ylab<-paste('E@-', sep='') else ylab<-paste(ylab, '@-', sep='')
		dmTrack <- DataTrack(exprgr[seqnames(exprgr)==prange[1] & (strand(exprgr)=='-' | strand(exprgr)=='*') & start(exprgr) >= as.numeric(prange[2]) & end(exprgr) <= as.numeric(prange[3])], name = ylab, background.title = "darkblue", type = c('p'),   legend = TRUE)
		trackList<-append(trackList, dmTrack)	
	}

	# start plotting
	# initial grapic dev.
	
	if(tofile){
		graphics.off()
		if(eps){
			setEPS()
			postscript(paste(main, '.', plotstrand,'.eps', sep=''), paper='special', width=25, height=10, horizontal=F, fonts=c("sans"))
		} else {
			pdf(paste(main, '.', plotstrand, '.pdf', sep=''), width=25, height=10)
		}
		# append general plotting params.
		# min.distance = 0, min.width = 0, to distinguish adjacent feature within the same group, not sure if there will be unexpected side effects or not.
		params <- append(params, list(fontsize=20, cex.title=1.5, main=main, min.distance = 0, min.width = 0)) 
		do.call(plotTracks,c(list(trackList), params)) 
		dev.off()
	} else {
		params <- append(params, list(main=main, min.distance = 0, min.width = 0)) 
		do.call(plotTracks,c(list(trackList), params)) 
	}
}

