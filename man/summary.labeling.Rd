\name{summary.labeling}
\alias{summary.labeling}

\title{
Summarize the results contained in an object of class \code{labeling}
}

\description{
This function summarizes the results contained in an object of the class \code{labeling}, giving both the estimated isotopic X abundances and the related standard errors.
}

\usage{
summary(fitted_abundances)
}


\arguments{
  \item{fitted_abundances}{
An object of class \code{labeling}. See \code{\link{main_labeling}} or \code{\link{find_abundance}} for the details.
}

}


\value{
A matrix summarizing the results contained in the object of class \code{labeling}. 
It has two rows: in the first are the estimated isotopic abundances of the labeling X isotope (\eqn{^2}H or \eqn{^{13}}C), in the second are the related standard errors coming from the \code{nls} fitting procedure. One column for each sample.
}