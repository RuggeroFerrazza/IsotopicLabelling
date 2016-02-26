\name{find_abundance}
\alias{find_abundance}

\title{
Fitting experimental isotopic patterns
}

\description{
Function that takes each of the provided experimental MS isotopic patterns, and fits the best theoretical pattern that reproduces it by iteratively changing the isotopic abundance of X (the labeling element).
}

\usage{
find_abundance(patterns, info, initial_abundance)
}

\arguments{
  \item{patterns}{
A matrix containing the experimental isotopic patterns in its columns (one column for each sample), with the first two representing the mass and the retention time of the corresponding peaks.
}
  \item{info}{
Named list, output of the \code{\link{isotopic_information}} function, containing useful information about the compound of interest and the isotopic distributions of its elements.
}
  \item{initial_abundance}{
Either NA, or a numeric vector of the same length as the number of samples, containing the initial guesses on the percentage isotopic abundance of the labeling isotope (denoted as X, it can be either \eqn{^{2}}H or \eqn{^{13}}C). If provided, numbers between 0 and 100. 
}
}

\value{
An object of class \code{\link{labeling}}, which is a list containing the results of the fitting procedure:
  \item{compound}{Character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution.}
  \item{best_estimate}{Numeric vector of length equal to the number of samples, containing the estimated percentage abundances of the labeling isotope X (either \eqn{^2}H or \eqn{^{13}}C). Numbers between 0 and 100.}
  \item{std_error}{Numeric vector with the standard errors of the estimated percentage abundances, output of the \code{nls} fitting procedure.}
  \item{dev_percent}{Numeric vector with the percentage deviations of the best fitted theoretical patterns to the provided experimental patterns.}
  \item{x_scale}{Numeric vector containing the m/z values relative to the signals of the experimental patterns.}
  \item{y_exp}{Matrix of normalized experimental isotopic patterns (one column for each sample). The highest signal of each pattern is set to 100.}
  \item{y_theor}{Matrix of normalized fitted theoretical isotopic patterns (one column for each sample). The highest signal of each pattern is set to 100.}
  \item{residuals}{Matrix of residuals: each column is the difference between (normalized) experimental and best fitted theoretical patterns.}
  \item{warnings}{Character vector of possible warnings from the fitting procedure.}
  }

\seealso{
\code{\link{isotopic_information}}
}