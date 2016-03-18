\name{plot.labelling}
\alias{plot.labelling}

\title{
Plot method for class "labelling" objects
}

\description{
It is a method for the plot generic. Depending on the parameter "type", three operating conditions are possible: by default (type="patterns") this function plots, for each sample in the class \code{labelling} object, its normalised experimental pattern superimposed to its fitted theoretical pattern.

By setting type="residuals", the function plots the residuals (the differences between experimental and best fitted theoretical patterns).

Finally, with type="summary" a summary plot is produced, showing the estimated percentage abundances with related standard errors.

If the parameter "saveplots" is set to TRUE, the plots are saved as a *.pdf file in the working directory.
}

\usage{
plot(fitted_abundances, type, saveplots)
}


\arguments{
  \item{fitted_abundances}{
Object of class \code{labelling}.}
  \item{type}{
The type of plots to produce: the default value is "patterns", but other two possible values are "residuals" and "summary" (see description above).
}
  \item{saveplots}{
A logical, with default value FALSE. If TRUE, the plots are saved to a single *.pdf file in the working directory.
}
}

\value{
One or more plots.
}