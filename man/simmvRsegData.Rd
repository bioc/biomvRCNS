\name{simmvRsegData}
\alias{simmvRsegData}
\title{
	Simulate exemplary segmentation data.
}
\description{
	Simulate exemplary segmentation data.
}
\usage{
simmvRsegData(nr = 50, nc = 10, N = c(3, 4), family = "norm", comVar = TRUE, segDisp = FALSE)
}

\arguments{
  \item{nr}{
	number of data points rows, representing feature number
}
  \item{nc}{
	number of columns, representing sample number
}
  \item{N}{
	number of change point per group
}
  \item{family}{
	family of distribution
}
  \item{comVar}{
	TRUE if a common covariance across samples for 'norm' data
}
  \item{segDisp}{
	TRUE if assuming a segment-wise dispersion parameter rather than a overall one
}
}

\value{
	a list object containing the start and end indices for new segments
	\item{x}{a nr*nc matrix of simulated data}
	\item{segStart}{a list of nc vectors of the starting indices of each segments corresponding to each x column}
	\item{segMean}{a list of nc vectors of each segment mean corresponding to each x column}
	\item{group}{the group affiliation of each x column}
}


\examples{
	x<-simmvRsegData()
}
