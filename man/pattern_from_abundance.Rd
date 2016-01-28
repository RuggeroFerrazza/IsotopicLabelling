\name{pattern_from_abundance}
\alias{pattern_from_abundance}

\title{
Computing theoretical isotopic patterns
}

\description{
Function that computes the theoretical pattern of the specified compound, given the abundance of the heaviest isotope of element X, contained in the compound. The function returns a vector of normalized theoretical intensities, where the maximum signal is set to 100; the related masses are contained in the info$target vector (from the input list).
}

\usage{
pattern_from_abundance(abundance, info)
}

\arguments{
  \item{abundance}{
Isotopic abundance of the heaviest isotope of the X element (either \eqn{^2}H or \eqn{^{13}}C). Number between 0 and 1.
}
  \item{info}{
A list with the required information, output of the \code{\link{isotopic_information}} function. 
}
}

\value{
A vector containing the normalized intensities of the theoretical pattern expected if X had the isotopic distribution specified at input.
}

\references{
The function makes use of the \code{\link[ecipex]{ecipex}} R package. 
}

\seealso{
\code{\link{isotopic_information}}, \code{\link[ecipex]{ecipex}}
}