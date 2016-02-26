\name{plot.labeling}
\alias{plot.labeling}

\title{
Plot method for labeling objects
}

\description{
It is a method for the plot generic. Depending on the parameter "type", three operating conditions are possible: by default (type="patterns") this function plots, for each sample contained in the object of class labeling, both its normalized experimental pattern and the best fitted theoretical pattern (maximum signals set to 100). The estimated percentage abundance of the labeling isotope (\eqn{^2}H or \eqn{^{13}}C) is also given. 

The function can also plot the residuals (type="residuals"), which are the differences between experimental and best fitted theoretical patterns.

Finally, a summary plot can be produced (type="summary") showing the estimated percentage abundances with related standard errors.

The "saveplots" parameter allows to save the plots as a single *.pdf document. 
}

\usage{
plot(fitted_abundances, type, saveplots)
}


\arguments{
  \item{fitted_abundances}{
Object of class \code{labeling}.}
  \item{type}{
The type of plots to produce: the default value is "patterns", which produces a series of plots (one for each sample) where the normalized experimental patterns are superimposed with their best fitted theoretical pattern.
}
  \item{saveplots}{
A logical, with default value FALSE. If TRUE, the plots are saved to a single *.pdf file in the working directory.
}
}

\value{
One or more plots.
}