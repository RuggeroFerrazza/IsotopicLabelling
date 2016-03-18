\name{IsotopicLabelling-package}
\alias{IsotopicLabelling-package}
\alias{IsotopicLabelling}
\docType{package}
\title{
\packageTitle{IsotopicLabelling}
}
\description{
The \code{IsotopicLabelling} package allows to analyse the isotopic patterns in MS data obtained following isotopic labelling experiments; by looking at the isotopic patterns, the package estimates the isotopic abundance of the stable isotope employed in the labelling experiment (either \eqn{^{2}}H or \eqn{^{13}}C) within the specified compound. 
}
\details{
Given a data frame of LC-MS or GC-MS peak intensities or areas (one column for each sample to analyse), the \code{\link{IsotopicLabelling}} package first extracts the isotopic patterns of the specified compound, and then performs an isotopic pattern analysis to estimate the isotopic abundance of the labelling isotope. This is obtained through a fitting procedure by which the abundance is iteratively changed and the resulting expected patterns are compared to the experimental one; the final estimated isotopic abundance is the value for which the theoretical pattern best reproduces the experimental one. The theoretical patterns are computed using the \code{\link[ecipex]{ecipex}} R package.

The isotopic pattern analysis can be divided into the following steps:
\enumerate{
\item Starting from a class \code{xcmsSet} object (from the \code{xcms} R package), generate a data frame of peak signal intensities or areas, with each column corresponding to a sample. This step can be avoided if the data frame is already available (obtained by other means);
\item Extract from the data frame the experimental isotopic patterns of the specified compound (one pattern for each sample). In the chemical formula of the compound, the element whose abundance is unknown is called "X";
\item Normalise the patterns and find the theoretical patterns that best fit them, by iteratively changing the isotopic abundance of the labelling isotope; 
\item Summarize the results.
}

}

\author{
\packageAuthor{IsotopicLabelling}

Maintainer: \packageMaintainer{IsotopicLabelling}
}



\examples{
data(xcms_obj)
peak_table <- table_xcms(xcms_obj)
fitted_abundances <- main_labelling(peak_table, compound="X40H77NO8P", labelling="C", 
mass_shift=0.05, RT=285, RT_shift=20, chrom_width=7, initial_abundance=NA)
summary(fitted_abundances)
plot(fitted_abundances, type="patterns", saveplots=FALSE)
plot(fitted_abundances, type="residuals", saveplots=FALSE)
plot(fitted_abundances, type="summary", saveplots=FALSE)
save_labelling(fitted_abundances)

# Other possible lipid compounds include:
# [PC34:1 + H]+. compound="X42H83NO8P", RT=475, chrom_width=10
# [TAG47:3 + NH4]+ (a minor species). compound="X50H94NO6", RT=891, chrom_width=7
}
