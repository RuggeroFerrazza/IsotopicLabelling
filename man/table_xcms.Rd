\name{table_xcms}
\alias{table_xcms}

\title{
Create a data frame of peak intensities or areas from an xcmsSet object
}

\description{
The function generates a data frame of MS peak intensities or areas starting from an object of class \code{xcmsSet} (from the \code{xcms} R package).
}

\usage{
table_xcms(xcms_obj)
}

\arguments{
  \item{xcms_obj}{
An object of class \code{xcmsSet}, output of the MS data processing with the \code{xcms} R package.
}
}

\value{
A data frame of experimental MS peak intensities or areas (one column for each sample), with its first two columns containing the \emph{m/z} and the retention times of the peaks.
}


\note{
The output data frame, required by other functions of the \code{\link{IsotopicLabelling}} R package, can be obtained in a number of other independent ways, such as through proprietary software of the vendor of the MS instrument. 
}


\seealso{
The \code{xcms} R package.
}

\examples{
data(xcms_obj)
peak_table <- table_xcms(xcms_obj)
}
