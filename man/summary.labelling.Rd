\name{summary.labelling}
\alias{summary.labelling}

\title{
Summarize the results contained in a class "labelling" object
}

\description{
This function summarizes the results contained in an object of class \code{labelling}, by showing the estimated percentage isotopic abundances of the labelling isotope, together with the related standard errors.
}

\usage{
summary(fitted_abundances)
}


\arguments{
  \item{fitted_abundances}{
An object of class \code{labelling}, output of either \code{\link{main_labelling}} or \code{\link{find_abundance}} functions.
}

}


\value{
A matrix summarizing the results contained in the class \code{labelling} object.
It has two rows: in the first one are the estimated percentage isotopic abundances of the labelling isotope X (\eqn{^2}H or \eqn{^{13}}C), in the second one are the related standard errors.
}