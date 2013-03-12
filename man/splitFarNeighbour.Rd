\name{splitFarNeighbour}
\alias{splitFarNeighbour}
\title{
	Split segments if long gaps exist between feature positions
}
\description{
	Split segments if long gaps exist between feature positions, due to low coverage or resolution.
}
\usage{
splitFarNeighbour(intStart = NULL, intEnd = NULL, xPos = NULL, xRange = NULL, maxgap = Inf, minrun = 1)
}

\arguments{
  \item{intStart}{
	indices of start for each segment 
}
  \item{intEnd}{
	indices of end for each segment 
}
  \item{xPos}{
	position vector, the distance of neighbouring features will be counted as point to point
}
  \item{xRange}{
	IRanges / GRanges object for the positions, the the distance of neighbouring features will be counted as end to start.
}
  \item{maxgap}{
	maximum distance between neighbouring features
}
  \item{minrun}{
	when splitting, the minimum length of the features spanning, which half will be ignored if shorter.
}
}

\value{
	a list object containing the start and end indices for new segments
	\item{IS}{the start indices for new segments}
	\item{IE}{the end indices for new segments}
}

\author{
Yang Du
}

\examples{
	set.seed(123)
	pos<-cumsum(rnbinom(20, size=10, prob=0.01))
	splitFarNeighbour(intStart=c(1, 10), intEnd=c(6, 18), xPos=pos, maxgap=1000)
}
