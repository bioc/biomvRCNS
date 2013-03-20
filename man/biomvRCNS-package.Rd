\name{biomvRCNS-package}
\alias{biomvRCNS-package}
\alias{biomvRCNS}
\docType{package}
\title{
biomvRCNS - Copy Number study and Segmentation models for multivariate biological data
}
\description{
In this package, a Hidden Semi Markov Model and one homogeneous segmentation based model are designed and implemented for genomic data, with the aim of assisting in transcripts detection using high throughput technology like RNA-seq or tiling array, and copy number analysis using array based CGH or targeted sequencing.
}
\details{
\tabular{ll}{
Package: \tab biomvRCNS\cr
Type: \tab Package\cr
Version: \tab 0.99.3\cr
Date: \tab 2013-03-20\cr
License: \tab GPL (>= 2)\cr
Imports: \tab methods, mvtnorm\cr
Depends: \tab IRanges, GenomicRanges, Gviz\cr
Suggests: \tab parallel, GenomicFeatures, dynamicTreeCut, Rsamtools, 	TxDb.Hsapiens.UCSC.hg19.knownGene\cr
}
}
\author{
Yang Du

Maintainer: Yang Du <duyang@fbn-dummerstorf.de>
}
\references{
	Guedon, Y. (2003). Estimating hidden semi-Markov chains from discrete sequences. Journal of Computational and Graphical Statistics, 12(3), 604-639.
	
	Piegorsch, W. W. (1990). Maximum likelihood estimation for the negative binomial dispersion parameter. Biometrics, 863-867.	

	Picard,F. et al. (2005) A statistical approach for array CGH data analysis. BMC Bioinformatics, 6, 27.
	
	Huber,W. et al. (2006) Transcript mapping with high density oligonucleotide tiling arrays. Bioinformatics, 22, 1963-1970. .

	Zhang, N. R. and Siegmund, D. O. (2007). A Modified Bayes Information Criterion with Applications to the Analysis of Comparative Genomic Hybridization Data. Biometrics 63 22-32.
	
	Robinson MD and Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data. Biostatistics, 9, 321-332
}

\keyword{ package }
\seealso{
  \code{\link{biomvRhsmm}}, 
  \code{\link{biomvRseg}}
}
