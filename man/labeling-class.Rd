\name{labeling-class}
\alias{labeling}

\title{
Class labeling
}

\description{
An object of the class labeling is a list containing the results of the analysis of MS data relative to isotopic labeling experiments.

These objects can be created by the functions \code{\link{find_abundance}} and \code{\link{main_labeling}}.
}

\details{
The methods available for class "labeling" objects are:
\describe{
	\item{summary}{to produce a summary of the results contained in the class labeling object.}
	\item{plot}{to create plots of the fitted isotopic patterns in the class labeling object.}
}

}

\seealso{
\code{\link{find_abundance}}, \code{\link{main_labeling}}, \code{\link{summary}}, \code{\link{plot}}
}