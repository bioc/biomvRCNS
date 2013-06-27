\name{biomvRGviz}
\alias{biomvRGviz}
\title{
	Plot segmentation result using Gviz
}
\description{
  This function could be called to plot segmentation output, together with the input signal and optional annotation. By default resulting image will be printed to file.
  The plot method for class \code{\link{biomvRCNS-class}} also calls this method. See the vignette for a more complete example.
}
\usage{
biomvRGviz(exprgr, gmgr = NULL, prange = NULL, regionID = NULL, seggr = NULL, plotstrand = "+", eps = TRUE, tofile = TRUE, ...)
}

\arguments{
  \item{exprgr}{ 
	a GRanges object with one numeric column for the segmentation input signal in its meta DataFrame
}
  \item{gmgr}{
	an optional GRanges object for the annotation, which must have one column named 'TYPE' in its meta DataFrame
}
  \item{prange}{
	an optional vector defining the scope of the plot, the first 3 elements must be formatted as c('seqname', 'start_position', 'end_position')
}
  \item{regionID}{
	a character for the name of the plotted region or gene name or other identifier, will be used in the title of the plot and the output file name
}
  \item{seggr}{
	a GRanges object for the segmentation output, which must have one column named 'STATE' in its meta DataFrame
}
  \item{plotstrand}{
	select which strand to plot, possible values are '+', '-', '*'
}
  \item{eps}{
	whether to output EPS file using postscript, if FALSE then PDF files for each sequence will be generated to the current working folder.
}
  \item{tofile}{
	whether to output graphics file, if FALSE then will plot on the active device.
}
  \item{...}{
	other arguments for \code{plot}, like \code{main} or \code{ylab}
}
}
\details{
	See the vignette for more details and examples.
}
\value{
	Plot graph on the active device or output to EPS/PDF file.
}

\author{
Yang Du
}


\examples{
	data(coriell)
	x<-coriell[coriell[,2]==1,]
	xgr<-GRanges(seqnames=paste('chr', x[,2], sep=''), IRanges(start=x[,3], width=1, names=x[,1]))
	values(xgr)<-DataFrame(x[,4:5], row.names=NULL)
	xgr<-xgr[order(xgr)]

	J<-2; maxk<-50
	# a uniform inital sojourn, not utilizing positional information
	soj<-list(J=J, maxk=maxk, type='gamma', d=cbind(dunif(1:maxk, 1, maxk), dunif(1:maxk, 1, maxk)))
	soj$D <- sapply(1:J, function(j) rev(cumsum(rev(soj$d[1:maxk,j]))))
	sample<-colnames(coriell)[4]
	runout<-hsmmRun(matrix(values(xgr)[,sample]), sample, xgr, soj, emis.type='norm')
	biomvRGviz(exprgr=xgr, seggr=runout$res, tofile=FALSE) 
	
}
