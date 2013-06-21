\name{biomvRCNS-class}
\Rdversion{1.1}
\docType{class}
\alias{biomvRCNS-class}
\alias{plot,biomvRCNS,ANY-method}
\alias{show,biomvRCNS-method}

\title{Class \code{"biomvRCNS"}}
\description{
	The default object class returned by \code{\link{biomvRhsmm}}, \code{\link{biomvRseg}} and \code{\link{biomvRmgmr}}
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("biomvRCNS", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"GRanges"}, with range information either from real positional data or just indices, with input data matrix stored in the meta columns.
    		 Additional meta columns for the estimated states and associated probabilities for each sample will also be appended following the input data matrix when using \code{biomvRhsmm}.} 
    \item{\code{res}:}{Object of class \code{"GRanges"} , each range represent one continuous segment identified, with sample name slot 'SAMPLE' and segment mean slot 'MEAN' stored in the meta columns  }
    \item{\code{param}:}{Object of class \code{"list"}, list of all parameters used in the corresponding model. }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "biomvRCNS", y = "ANY")}: ... }
    \item{show}{\code{signature(object = "biomvRCNS")}: ... }
	 }
}
\examples{
showClass("biomvRCNS")
}
\keyword{classes}
