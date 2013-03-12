\name{biomvRseg}
\alias{biomvRseg}
\title{
	Homogeneous segmentation of multi-sample genomic data
}
\description{
 The function will perform a two stage segmentation on multi-sample genomic data from array experiment or high throughput sequencing data.
}
\usage{
biomvRseg(x, maxk=NULL, maxbp=NULL, maxseg=NULL, xPos=NULL, xRange=NULL, usePos='start', family='norm', penalty='BIC', twoStep=TRUE, segDisp=FALSE, useMC=FALSE, useSum=TRUE, comVar=TRUE, maxgap=Inf, tol=1e-06, grp=NULL, clusterm=NULL, na.rm=TRUE)
}

\arguments{
  \item{x}{
	input data matrix, or a \code{GRanges} object with input stored in the meta \code{DataFrame}
}
  \item{maxk}{
	maximum length of a segment
}
  \item{maxbp}{
	maximum length of a segment in bp, given positional information specified in \code{xPos} / \code{xRange} / or \code{x}
}
  \item{maxseg}{
	maximum number of segment the function will try
}

 \item{xPos}{
	a vector of positions for each \code{x} row
}
  \item{xRange}{
	a \code{IRanges}/\code{GRanges} object, same length as \code{x} rows
}
  \item{usePos}{
	character value to indicate whether the 'start', 'end' or 'mid' point position should be used
}
  \item{family}{
	family of \code{x} distribution, only the following types are supported: 'norm', 'nbinom', 'pois'
}
  \item{penalty}{
  	 penalty method used for determining the optimal number of segment using likelihood, possible values are 'none','AIC','AICc','BIC','SIC','HQIC', 'mBIC'
}
  \item{twoStep}{
	TRUE if a second stage merging will be performed after the initial group segmentation
}
  \item{segDisp}{
	TRUE if a segment-wise estimation of dispersion parameter rather than using a overall estimation
}
  \item{useMC}{
	TRUE if \code{\link{mclapply}} should be used to speed up the calculation for nbinom dispersion estimation
}
  \item{useSum}{
	TRUE if using grand sum across sample / x columns, like in the \code{tilingArray} solution
}
  \item{comVar}{
	TRUE if assuming common variance across samples (x columns)
}
  \item{maxgap}{
	max distance between neighbouring feature to consider a split
}
  \item{tol}{
	tolerance level of the likelihood change to determining the termination of the EM run
}
  \item{grp}{
	vector of group assignment for each sample, with a length the same as columns in the data matrix, samples within each group would be processed simultaneously if a multivariate emission distribution is available
}
  \item{clusterm}{
	clustering method for prior grouping, possible values are 'ward','single','complete','average','mcquitty','median','centroid'
}
  \item{na.rm}{
	\code{TRUE} if \code{NA} value should be ignored 
}
}
\details{
	A homogeneous segmentation algorithm, using dynamic programming like in \code{tilingArray}; however capable of handling count data from sequencing.
}
\value{
	A \code{\link{biomvRCNS}} object:
    \item{\code{x}:}{Object of class \code{"GRanges"}, with range information either from real positional data or just indices, with input data matrix stored in the meta columns.} 
    \item{\code{res}:}{Object of class \code{"GRanges"} , each range represent one continuous segment identified, with sample name slot 'SAMPLE' and segment mean slot 'MEAN' stored in the meta columns  }
    \item{\code{param}:}{Object of class \code{"list"}, list of all parameters used in the model run. }
}
\references{
	Piegorsch, W. W. (1990). Maximum likelihood estimation for the negative binomial dispersion parameter. Biometrics, 863-867.	

	Picard,F. et al. (2005) A statistical approach for array CGH data analysis. BMC Bioinformatics, 6, 27.
	
	Huber,W. et al. (2006) Transcript mapping with high density oligonucleotide tiling arrays. Bioinformatics, 22, 1963-1970. .

	Zhang, N. R. and Siegmund, D. O. (2007). A Modified Bayes Information Criterion with Applications to the Analysis of Comparative Genomic Hybridization Data. Biometrics 63 22-32.

	Robinson MD and Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data. Biostatistics, 9, 321-332

}

\seealso{
\code{\link{biomvRhsmm}}
}
\examples{
	data(coriell)
	xgr<-GRanges(seqnames=paste('chr', coriell[,2], sep=''), IRanges(start=coriell[,3], width=1, names=coriell[,1]))
	values(xgr)<-DataFrame(coriell[,4:5], row.names=NULL)
	xgr<-xgr[order(xgr)]
	resseg<-biomvRseg(x=xgr, maxbp=4E4, maxseg=10, family='norm', grp=c(1,2))
}
