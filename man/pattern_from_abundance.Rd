\name{pattern_from_abundance}
\alias{pattern_from_abundance}

\title{
Compute theoretical isotopic patterns
}

\description{
Function that computes the theoretical pattern of the specified compound, given the abundance of the labeling isotope X. The function returns a vector of normalised theoretical intensities, where the maximum signal is set to 100; the related masses are in the "target" vector contained in the input list.
}

\usage{
pattern_from_abundance(abundance, info)
}

\arguments{
  \item{abundance}{
Isotopic abundance of the labeling isotope X (either \eqn{^2}H or \eqn{^{13}}C); number between 0 and 1.
}
  \item{info}{
Named list containing isotopic information, output of the \code{\link{isotopic_information}} function.
}
}

\value{
A vector representing the normalised isotopic pattern of the compound of interest, corresponding to the specified isotopic distribution.
}

\references{
The function makes use of the \code{\link[ecipex]{ecipex}} R package. 
}

\seealso{
\code{\link{isotopic_information}}, \code{\link[ecipex]{ecipex}}
}