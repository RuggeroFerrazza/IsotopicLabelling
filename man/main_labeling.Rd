\name{main_labeling}
\alias{main_labeling}

\title{
Estimating isotopic abundances from experimental MS data
}

\description{
The main function of the \code{\link{IsotopicLabeling}} package, it estimates the isotopic abundance of the labeling element in the specified compound, given a table of MS peak intensities or areas.
}

\usage{
main_labeling(peak_table, compound, labeling, mass_shift, 

 RT, RT_shift, chrom_width, initial_abundance)
}


\arguments{
  \item{peak_table}{
Data frame containing experimental peak intensities or areas relative to a set of samples (one for each column), with the first two columns representing the m/z and the retention time of the peaks.

If not already available, the table can be obtained from an \code{xcmsSet} object (from the \code{xcms} R package) through the \code{\link{table_xcms}} function.
}
  \item{compound}{
Character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted). 
}
  \item{labeling}{
Character, either "H" or "C", specifying which is the labeling element. 
}
  \item{mass_shift}{
Maximum mass shift between theoretical and experimental mass. In other words, the mass accuracy.
}
  \item{RT}{
Expected retention time of the compound of interest (in sec).
}
  \item{RT_shift}{
Maximum time shift between expected and experimental retention time of the peaks (in sec).
}
  \item{chrom_width}{
An estimate of the peak chromatographic width (in sec). 
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
  \item{y_exp}{Matrix of normalised experimental isotopic patterns (one column for each sample). The highest signal of each pattern is set to 100.}
  \item{y_theor}{Matrix of normalised fitted theoretical isotopic patterns (one column for each sample). The highest signal of each pattern is set to 100.}
  \item{residuals}{Matrix of residuals: each column is the difference between (normalised) experimental and best fitted theoretical patterns.}
  \item{warnings}{Character vector of possible warnings from the fitting procedure.}
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