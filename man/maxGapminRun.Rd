\name{maxGapminRun}
\alias{maxGapminRun}
\title{
	Max-gap-min-run algorithm for 2 states segmentation
}
\description{
	A custom Max-gap-min-run implementation using physical position for gap and run length calculation.
}
\usage{
maxGapminRun(x, xPos = NULL, xRange = NULL, cutoff = NULL, q = 0.9, high=TRUE, minrun = 5, maxgap = 2, splitLen = Inf, na.rm=TRUE)
}

\arguments{
  \item{x}{
	a numeric vector for the input signal
}
  \item{xPos}{
	a numeric vector, same length as x, carrying positional information for each element of x 
}
  \item{xRange}{
	an \code{IRanges} object,  same length as x, carrying range information for each element of x 
}
  \item{cutoff}{
	numeric value used as cut-off, optional if \code{q} is specified
}
  \item{q}{
	numeric value used to derive cut-off of x, as the \code{q} quantile of x , optional if \code{cutoff} is specified
}
  \item{high}{
	TRUE if the \code{cutoff} or \code{q} here is the lower bound and values greater than the threshold are considered 
}
  \item{minrun}{
  	minimum run length for the resulting segments
}
  \item{maxgap}{
	numeric value, maximum gap length between regions with TRUE signal.
}
  \item{splitLen}{
	numeric value, maximum length of segments, split if too long
}
  \item{na.rm}{
	\code{TRUE} if \code{NA} value should be ignored 
}
}
\details{
	A custom Max-gap-min-run implementation using physical position for gap and run length calculation.
}
\value{
	a list of segment starts and ends indices
  \item{IS}{the start index for each segment}
  \item{IE}{the end index for each segment}
  \item{CUTOFF}{the cutoff value used in the run}
  \item{MG}{the parameter value for \code{maxgap}}  
  \item{MR}{the parameter value for \code{minrun}}
  \item{SPL}{the parameter value for \code{splitLen}}
}

\author{
Yang Du
}

\seealso{
\code{\link{biomvRhsmm}}
\code{\link{biomvRseg}}
\code{\link{biomvRmgmr}}
}

\examples{
	x<-rpois(50, 10)
	xpos<-rnorm(50, 300, 100)
	xpos<-xpos[order(xpos)]
	maxGapminRun(x, xpos, cutoff=9.5, maxgap=30, minrun=100)
}

