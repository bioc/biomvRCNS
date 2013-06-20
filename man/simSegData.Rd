\name{simSegData}
\alias{simSegData}
\title{
	Simulate exemplary segmentation data.
}
\description{
	Simulate exemplary segmentation data.
}
\usage{
simSegData(nseg=10, J=3, soj, emis, seed=1234, toPlot=FALSE)
}

\arguments{
  \item{nseg}{
	size of initial segments pool
}
  \item{J}{
	states number
}
  \item{soj}{
	a list object containing sojourn settings
}
  \item{emis}{
	a list object containing emission settings
}
  \item{seed}{
	seed for simulation
}
  \item{toPlot}{
	whether to output a pdf image of the simulated series
}
}

\value{
	a list object containing the simulated data and the segment info
	\item{E}{a numeric vector of the simulated data serie}
	\item{L}{a vector of the length for each continuous segment}
	\item{S}{a vector of state assignment for each segment}
	\item{pdf}{the name of the output pdf file if any}
}


\examples{
	soj<-list(type='pois', lambda=c(200, 100, 10))
	emis<-list(type='pois', lambda=1:3)
	simSegData(soj=soj, emis=emis)
}
