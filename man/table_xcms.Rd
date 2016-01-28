\name{table_xcms}
\alias{table_xcms}

\title{
Create a table of peak intensities or areas from an xcmsSet object
}

\description{
The function generates a table of peak intensities or areas from an object of the class xcmsSet (from the \code{xcms} R package).
}

\usage{
table_xcms(xcms_obj)
}

\arguments{
  \item{xcms_obj}{
An object of the xcmsSet class, output of the data processing with the \code{xcms} R package.
}
}

\value{
A data frame of experimental peak intensities or areas (one column for each sample), with its first two columns containing the m/z and the retention times of the peaks.
}


\note{
The output table, required by the function \code{\link{isotopic_pattern}} or \code{\link{main_labeling}}, can be obtained in a number of other independent ways, such as through proprietary software of the vendor of the MS instrument. 
}


\seealso{
The \code{xcms} R package.
}

\examples{
data(xcms_obj)
peak_table <- table_xcms(xcms_obj)
}
