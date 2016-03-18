\name{main_labelling}
\alias{main_labelling}

\title{
Estimate isotopic abundances from experimental MS isotopic patterns
}

\description{
The main function of the \code{\link{IsotopicLabelling}} package, it estimates the isotopic abundance of the labelling isotope in the specified compound, given a data frame containing MS peak intensities or areas.
}

\usage{
main_labelling(peak_table, compound, labelling, mass_shift, 

 RT, RT_shift, chrom_width, initial_abundance)
}


\arguments{
  \item{peak_table}{
Data frame of experimental MS peak intensities or areas (one column for each sample), with the first two columns representing \emph{m/z} and retention times of the peaks.

If not already available, the table can be obtained from an \code{xcmsSet} object (output of the \code{xcms} R package) through the \code{\link{table_xcms}} function.
}
  \item{compound}{
Character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted). 
}
  \item{labelling}{
Character, either "H" or "C", specifying the labelling element. 
}
  \item{mass_shift}{
Maximum difference between theoretical and experimental mass. In other words, the mass accuracy.
}
  \item{RT}{
Expected retention time of the compound of interest.
}
  \item{RT_shift}{
Maximum difference between expected and experimental retention time of the peaks.
}
  \item{chrom_width}{
An estimate of the chromatographic peak width.
}
  \item{initial_abundance}{
Either NA, or a numeric vector of length equal to the number of samples, with the initial guesses on the percentage isotopic abundance of the labelling isotope (denoted as X, it can be either \eqn{^{2}}H or \eqn{^{13}}C). If provided, numbers between 0 and 100. 
}
}

\value{
An object of class \code{\link{labelling}}, which is a list containing the results of the fitting procedure:
  \item{compound}{Character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted).}
  \item{best_estimate}{Numeric vector of length equal to the number of samples, containing the estimated percentage abundances of the labelling isotope X (either \eqn{^2}H or \eqn{^{13}}C). Numbers between 0 and 100.}
  \item{std_error}{Numeric vector with the standard errors of the estimates, output of the \code{nls} fitting procedure.}
  \item{dev_percent}{Numeric vector with the percentage deviations between best fitted and related experimental patterns.}
  \item{x_scale}{Numeric vector containing the m/z values relative to the signals of the experimental patterns.}
  \item{y_exp}{Matrix of normalised experimental isotopic patterns (one column for each sample). The most intense signal of each pattern is set to 100.}
  \item{y_theor}{Matrix of normalised fitted theoretical isotopic patterns (one column for each sample). The most intense signal of each pattern is set to 100.}
  \item{residuals}{Matrix of residuals: each column is the difference between experimental and best fitted theoretical patterns.}
  \item{warnings}{Character vector with possible warnings from the \code{nls} fitting procedure.}
  }


\note{
This function is made up of three important functions of the \code{\link{IsotopicLabelling}} package which are, in order: \code{\link{isotopic_information}}, \code{\link{isotopic_pattern}} and \code{\link{find_abundance}}. The user can therefore choose whether to carry out all the processing in a single step using this compact function, or to separately perform the three successive steps without using \code{\link{main_labelling}}.
}

\seealso{
\code{\link{table_xcms}}
}

\examples{
data(xcms_obj)
peak_table <- table_xcms(xcms_obj)
fitted_abundances <- main_labelling(peak_table, compound="X40H77NO8P", labelling="C", 
mass_shift=0.05, RT=285, RT_shift=20, chrom_width=7, initial_abundance=NA)

}