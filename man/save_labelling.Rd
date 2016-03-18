\name{save_labelling}
\alias{save_labelling}

\title{
Save class "labelling" objects to *.csv files
}

\description{
This function allows to save the results of the isotopic pattern analysis to a *.csv file.
}

\usage{
save_labelling(fitted_abundances, path)
}


\arguments{
  \item{fitted_abundances}{
Object of class \code{labelling}. It contains the results of the isotopic pattern analysis.
}
  \item{path}{
The directory where to save the *.csv file. If not specified, it is saved in the working directory. 
}

}


\value{
The "COMPOUND_Estimated_Abundances.csv" file, containing the results of the isotopic pattern analysis. For each sample (one for each row) there are four columns:
\enumerate{
\item The estimated percentage abundance of the labelling isotope (either \eqn{^2}H or \eqn{^{13}}C);
\item The related standard error coming from the \code{nls} fitting procedure;
\item The percentage deviation between theoretical and experimental isotopic patterns;
\item The outcome message from the fitting procedure, to undersand whether there have been any convergence problems.
}
}

\seealso{
\code{\link{main_labelling}}
}