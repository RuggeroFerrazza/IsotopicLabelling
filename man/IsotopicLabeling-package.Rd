\name{IsotopicLabeling-package}
\alias{IsotopicLabeling-package prova2}
\alias{IsotopicLabeling}
\docType{package}
\title{
\packageTitle{IsotopicLabeling}
}
\description{
\packageDescription{IsotopicLabeling}
}
\details{
Given a table of peak intensities or areas (one column for each sample), the \code{\link{IsotopicLabeling}} package first extracts the isotopic patterns of the specified compound, and then analyzes them in order to estimate the isotopic abundance of the labeling isotope. Such estimate is obtained through a fitting procedure by which the isotopic abundance of the labeling isotope is iteratively changed, and each time the resulting expected theoretical pattern is compared to the experimental one; the theoretical patterns are computed through the \code{\link[ecipex]{ecipex}} package. The final estimated isotopic abundance corresponds to the theoretical pattern that best reproduces the experimental pattern. 

The analysis performed by the \code{\link{IsotopicLabeling}} package can be divided into the following steps:
\enumerate{
\item Starting from an object of the xcmsSet class (from the \code{xcms} R package), create a table of peak signal intensities or areas, with each column corresponding to a sample. To this aim the \code{\link{table_xcms}} function is used, but this step can be skipped as well if the table is already available;
\item Extract from the table the experimental isotopic patterns of the specified compound (one pattern for each sample). In the chemical formula of the compound, the element whose abundance is unknown should be denoted as "X", such as "NX4" for the ammonium ion with all 4 hydrogen atoms having unknown abundance;
\item Normalize the experimental patterns and find the theoretical patterns that best fit them, by iteratively changing the isotopic abundance of element X; 
\item Provide the results. Both plots and summary tables can be created and saved.
}

}

\author{
\packageAuthor{IsotopicLabeling}

Maintainer: \packageMaintainer{IsotopicLabeling}
}



\examples{
data(xcms_obj)
peak_table <- table_xcms(xcms_obj)
fitted_abundances <- main_labeling(peak_table, compound="X40H77NO8P", labeling="C", 
mass_shift=0.05, RT=285, RT_shift=20, chrom_width=7, initial_abundance=NA)
summary(fitted_abundances)
plot(fitted_abundances, path=NA, what="patterns")
plot(fitted_abundances, path=NA, what="residuals")
save_labeling(fitted_abundances)

# Other possible lipid compounds include:
# [PC34:1 + H]+. compound="X42H83NO8P", RT=475, chrom_width=10
# [TAG47:3 + NH4]+ (a minor species). compound="X50H94NO6", RT=891, chrom_width=7
}
