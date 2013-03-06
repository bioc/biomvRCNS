\name{regionSegAlphaNB}
\alias{regionSegAlphaNB}
\title{
	Estimate matrix of dispersion parameter alpha (size) used in \code{\link{regionSegCost}} for negative binomial distributed \code{x}.
}
\description{
	Estimate matrix of dispersion parameter alpha (size) used in \code{\link{regionSegCost}} for negative binomial distributed \code{x}.
}
\usage{
regionSegAlphaNB(x, maxk = NULL, segs = NULL, useMC = FALSE, tol=1e-06)
}

\arguments{
  \item{x}{
	The input data matrix or vector
}
  \item{maxk}{
	Maximum number of index to search forward
}
  \item{segs}{
	Starting indices (excluding 1) for the candidate segments, for the second stage model, \code{maxk} will be overridden with \code{length(segs)+1}.
}
  \item{useMC}{
	TRUE if \code{\link{mclapply}} should be used to speed up the calculation
}
  \item{tol}{
	tolerance level for the convergence criteria in the maximum likelihood estimation of negative binomial distribution dispersion parameter.
}
}
\details{
	Estimate matrix of dispersion parameter alpha (size) used in \code{\link{regionSegCost}} for negative binomial distributed \code{x}.
}
\value{
	Matrix with \code{maxk} rows and \code{nrow(x)} columns, or a \code{length(segs)+1} square matrix for the second stage model.
}
\references{
	Piegorsch, W. W. (1990). Maximum likelihood estimation for the negative binomial dispersion parameter. Biometrics, 863-867.	
	
	Robinson MD and Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data. Biostatistics, 9, 321-332
}


\seealso{
\code{\link{regionSegCost}}
}

\examples{
	x<-matrix(rnbinom(120, size=0.05, mu=20), ncol=3)
	Aa<-regionSegAlphaNB(x, maxk=20)
	dim(Aa) # [1] 20 40
	Ab<-regionSegAlphaNB(x, segs=as.integer(c(3, 6, 12, 30)))
	dim(Ab) # [1] 5 5
}

