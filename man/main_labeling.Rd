\name{main_labeling}
\alias{main_labeling}

\title{
Estimating isotopic abundances from experimental LC-MS data
}

\description{
The main function of the \code{\link{IsotopicLabeling}} package, it estimates the isotopic abundance of the labeling element in the specified compound, given a table of peak intensities or areas as its input.
More samples can be contained in the table (one for each column).
}

\usage{
main_labeling(peak_table, compound, labeling, mass_shift, 

 RT, RT_shift, chrom_width, initial_abundance)
}


\arguments{
  \item{peak_table}{
Data frame containing experimental peak intensities or areas relative to a set of samples (one sample for each column), with the first two columns representing the m/z and the retention time of the peaks.

If not available, the table can be obtained from the function \code{\link{table_xcms}}, starting from an xcmsSet object (from the \code{xcms} R package).
}
  \item{compound}{
Character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted). 
}
  \item{labeling}{
Character, either "H" or "C", specifying which is the labeling element. 
}
  \item{mass_shift}{
Maximum shift allowed in the mass range between theoretical and experimental mass. In other words, the accuracy of the measurements.
}
  \item{RT}{
Expected retention time of the compound of interest (in sec).
}
  \item{RT_shift}{
Maximum shift allowed in the retention time range (in sec).
}
  \item{chrom_width}{
An estimate of the chromatographic width of the peaks (in sec). 
}
  \item{initial_abundance}{
Either NA, or a numeric vector of the same length as the number of samples, containing the initial guesses on the isotopic abundance of the labeling isotope (denoted as X, it can be either \eqn{^2}H or \eqn{^{13}}C). If provided, numbers between 0 and 1. 
}
}

\value{
An object of the class \code{\link{labeling}}, which is a list containing the results from the fitting procedure:
  \item{compound}{Character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution.}
  \item{Best_estimate}{Numeric vector of length equal to the number of samples, containing the estimated abundances of the labeling isotope (either \eqn{^2}H or \eqn{^{13}}C, and denoted as X). Numbers between 0 and 1.}
  \item{std_error}{Numeric vector with the standard errors of the estimated abundances; output of the \code{nls} fitting procedure.}
  \item{dev_percent}{Numeric vector with the percentage deviations of the best fitted theoretical patterns to the provided experimental patterns.}
  \item{x_scale}{Numeric vector containing the m/z values relative to the signals of the experimental patterns}
  \item{y_exp}{Matrix of normalized experimental isotopic patterns (one column for each sample); for each sample, the highest signal is set to 100.}
  \item{y_theor}{Matrix of normalized fitted theoretical isotopic patterns (one column for each sample); for each sample, the highest signal is set to 100.}
  \item{residuals}{Matrix containing the residuals: each column represents the difference between the observed normalized experimental pattern of a sample and its best fitted (normalized) theoretical pattern.}
  \item{warnings}{Character vector of possible warnings coming from the fitting procedure.}
  }


\note{
This function is built up of other functions of the \code{\link{IsotopicLabeling}} package which are, in order:   \code{\link{isotopic_information}}, \code{\link{isotopic_pattern}} and \code{\link{find_abundance}}. The user can therefore choose to carry out all the processing in a single step using this compact function. 
}

\seealso{
\code{\link{table_xcms}}
}

\examples{
data(xcms_obj)
peak_table <- table_xcms(xcms_obj)
fitted_abundances <- main_labeling(peak_table, compound="X40H77NO8P", labeling="C", 
mass_shift=0.05, RT=285, RT_shift=20, chrom_width=7, initial_abundance=NA)

}