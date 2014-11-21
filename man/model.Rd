\name{model}
\docType{data}
\alias{model}

\title{Components of the non-parametric mixture moodel fitted on Ago2 PAR-CLIP data}

\description{The non-parametric mixture model was fit on the entire Ago2 public available PAR-CLIP dataset (Kishore et al.) using the \code{fitMixtureModel} function.
}

\usage{
data(model)
model

}
\format{List of 5 items containing the estimated mixing coefficients and model densities. See the help page of the \link{fitMixtureModel} function for a detailed description of the output.
}

\references{
	Kishore et al. A quantitative analysis of CLIP methods for identifying binding sites of RNA-binding proteins. Nat Methods (2011) vol. 8 (7) pp. 559-64 \url{http://www.nature.com/nmeth/journal/v8/n7/full/nmeth.1608.html}
}

\seealso{
\code{\link{fitMixtureModel}}
}

\keyword{datasets}