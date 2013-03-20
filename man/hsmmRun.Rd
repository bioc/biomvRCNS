\name{hsmmRun}
\alias{hsmmRun}
\title{
 Estimating the most likely state sequence using Hidden Semi Markov Model
}
\description{
 This is the working horse of the \code{biomvRhsmm}
}
\usage{
hsmmRun(x, xid = "sampleid", xRange, soj, emis.type = "norm", q.alpha = 0.05, r.var = 0.75, cMethod='F-B', maxit = 1, maxgap = Inf, tol = 1e-06, avg.m='median', trim=0,na.rm = TRUE)
}

\arguments{
  \item{x}{
	input data matrix or vector, ordered wrt. position
}
  \item{xid}{
	name of the sample
}
  \item{xRange}{
	a \code{IRanges}/\code{GRanges} object, same length as \code{x} rows
}
  \item{soj}{
	a list object containing the relevant sojourn distribution parameters
}
  \item{emis.type}{
	type of the emission distribution, only the following types are supported: 'norm', 'mvnorm', 'pois', 'nbinom', 'mvt', 't'
}
  \item{q.alpha}{
	a quantile factor controlling the estimated prior for the mean of the emission of each states, \code{seq(from=q.alpha, to=1-q.alpha, length.out=J)}
}
  \item{r.var}{
	a ratio factor controlling the estimated prior for the variance / covariance structure of each states. 
	A value larger than 1 tend to allow larger variation in extreme states;  a value smaller than 1 will decrease the probability of having extreme states
}
  \item{cMethod}{
	C algorithm used for the most likely state sequence, 'F-B' or 'Viterbi'
}
 \item{maxit}{
	max iteration of the EM run with Forward-Backward algorithm
}
  \item{maxgap}{
	max distance between neighbouring feature to consider a split
}
  \item{tol}{
	tolerance level of the likelihood change to terminate the EM run
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
	The function fits a Hidden-semi Markov model for the input data matrix / vector, which should contains ordered data from a continuous region on one chromosome. 
	The model will start with flat prior for the initial state probability and transition probability, while emission parameter for each state will be estimated using different quantiles of the input controlled by argument \code{q.alpha} and \code{r.var}. 
	Argument for the sojourn density should be provided via the list object \code{soj}, which is either initialized as flat prior or estimated from other data in a previous call.
	The positional information in the \code{xRange} is used for the optional spiting of physically distant features and construction of returning \code{GRanges} object \code{res}.
}
\value{
	a list object, 
    \item{\code{yhat}}{the most likely state sequence, same length as \code{x} rows number}
    \item{\code{res}}{Object of class \code{"GRanges"} , each range represent one continuous segment identified, with sample name slot 'SAMPLE', estimated state slot 'STATE' and segment mean slot 'MEAN' stored in the meta columns}
}

\author{
Yang Du
}


\references{
	Guedon, Y. (2003). Estimating hidden semi-Markov chains from discrete sequences. Journal of Computational and Graphical Statistics, 12(3), 604-639.
}

\seealso{
\code{\link{biomvRhsmm}}
}

\examples{
	data(coriell)
	# select only chr1 
	x<-coriell[coriell[,2]==1,]
	xgr<-GRanges(seqnames=paste('chr', x[,2], sep=''), IRanges(start=x[,3], width=1, names=x[,1]))
	values(xgr)<-DataFrame(x[,4:5], row.names=NULL)
	xgr<-xgr[order(xgr)]

	J<-2 ; maxk<-50
	# a uniform initial sojourn, not utilizing positional information, just the index
	soj<-list(J=J, maxk=maxk, type='gamma', d=cbind(dunif(1:maxk, 1, maxk), dunif(1:maxk, 1, maxk)))
	soj$D <- sapply(1:J, function(j) rev(cumsum(rev(soj$d[1:maxk,j]))))
	# run 1 sample only, Coriell.13330
	sample<-colnames(coriell)[5]
	runout<-hsmmRun(matrix(values(xgr)[,sample]), sample, xgr, soj, emis.type='norm')
}

\keyword{ hsmm }
