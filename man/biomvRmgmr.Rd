\name{biomvRmgmr}
\alias{biomvRmgmr}
\title{
	Batch process multiple sequences and samples using max-gap-min-run algorithm for 2 states segmentation
}
\description{
	This is a wrapper function for batch processing multiple sequences and samples using max-gap-min-run algorithm for 2 states segmentation
}
\usage{
biomvRmgmr(x, xPos=NULL, xRange=NULL, usePos='start', cutoff=NULL, q=0.9, high=TRUE, minrun=5, maxgap=2, splitLen=Inf, poolGrp=FALSE, grp=NULL, cluster.m=NULL, avg.m='median', trim=0,na.rm=TRUE)
}

\arguments{
  \item{x}{
	input data matrix, or a \code{GRanges} object with input stored in the meta \code{DataFrame}, assume ordered.
}
 \item{xPos}{
	a vector of positions for each \code{x} row
}
  \item{xRange}{
	a \code{IRanges}/\code{GRanges} obejct, same length as \code{x} rows
}
  \item{usePos}{
	character value to indicate whether the 'start', 'end' or 'mid' point position should be used
}
  \item{cutoff}{
	threshold level above which is considered extreme
}
  \item{q}{
	relative quantile threshold level instead of absolute value for the cutoff
}
  \item{high}{
	TRUE if the \code{cutoff} or \code{q} here is the lower bound and values greater than the threshold are considered 
}
  \item{minrun}{
  	minimum run length for the resulting segments
}
  \item{maxgap}{
	maximum genomic distance below which two adjacent qualified tiles can be joined
}
  \item{splitLen}{
	numeric value, maximum length of segments, split if too long
}
  \item{poolGrp}{
	TRUE if samples within the same group should be pooled using median for each feature
}
  \item{grp}{
	vector of group assignment for each sample, with a length the same as columns in the data matrix, samples within each group would be processed simultaneously if a multivariate emission distribution is available
}
  \item{cluster.m}{
	clustering method for prior grouping, possible values are 'ward','single','complete','average','mcquitty','median','centroid'
}
  \item{avg.m}{
	method to calculate average value for each segment, 'median' or 'mean' possibly trimmed
}
  \item{trim}{
	the fraction (0 to 0.5) of observations to be trimmed from each end of x before the mean is computed. Values of trim outside that range are taken as the nearest endpoint.
}
  \item{na.rm}{
	\code{TRUE} if \code{NA} value should be ignored 
}
}
\details{
		This is the batch function to apply \code{\link{maxGapminRun}} multiple sequence.
}
\value{
	A \code{\link{biomvRCNS}} object:
    \item{\code{x}:}{Object of class \code{"GRanges"}, with range information either from real positional data or just indices, with input data matrix stored in the meta columns.} 
    \item{\code{res}:}{Object of class \code{"GRanges"} , each range represent one continuous segment identified, with sample name slot 'SAMPLE' and segment mean slot 'MEAN' stored in the meta columns  }
    \item{\code{param}:}{Object of class \code{"list"}, list of all parameters used in the model run. }
}

\author{
Yang Du
}

\seealso{
\code{\link{biomvRhsmm}}
\code{\link{maxGapminRun}}
}
\examples{
	data(coriell)
	xgr<-GRanges(seqnames=paste('chr', coriell[,2], sep=''), IRanges(start=coriell[,3], width=1, names=coriell[,1]))
	values(xgr)<-DataFrame(coriell[,4:5], row.names=NULL)
	xgr<-xgr[order(xgr)]
	resseg<-biomvRmgmr(x=xgr, minrun=3000, maxgap=1500, q=0.9, grp=c(1,2))
}
