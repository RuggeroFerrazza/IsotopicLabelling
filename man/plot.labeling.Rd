\name{plot.labeling}
\alias{plot.labeling}

\title{
Plot method for labeling objects
}

\description{
It is a method for the plot generic. Depending on the parameter "what", two operating conditions are possible: by default this function plots, for each sample contained in the object of class labeling, both its normalized experimental pattern and the best fitted theoretical pattern (maximum signals set to 100). The estimated abundance of the labeling isotope (\eqn{^2}H or \eqn{^{13}}C) is also given. A further plot is produced, showing the estimated abundances with related standard errors.

Alternatively, the function can plot the residuals (the differences between experimental and best fitted theoretical patterns).

A copy of the plots is saved to the specified directory as *.png files. 
}

\usage{
plot(fitted_abundances, path, what)
}


\arguments{
  \item{fitted_abundances}{
Object of class \code{labeling}.}
  \item{path}{
The directory where to save the plots. If not specified, a folder named after the chemical formula of the compound is created in the working directory, and the plots are saved there.
}
  \item{what}{
A character vector specifying what to plot, either "patterns" (the default value) or "residuals". In the first case, the normalized experimental patterns are plotted together with their best fitted theoretical patterns. In the second case, the residuals are plotted.
}
}

\value{
A series of plots. They are all printed to the screen, but a copy of them is also saved to the specified directory as *.png files.
}