\name{isotopic_pattern}
\alias{isotopic_pattern}

\title{
Extract experimental isotopic patterns from a table of MS peaks
}

\description{
Function that extracts the experimental isotopic patterns of a specified compound from a data frame containing MS peak intensities or areas.
}

\usage{
isotopic_pattern(peak_table, info, mass_shift, RT, RT_shift, chrom_width)
}


\arguments{
  \item{peak_table}{
Data frame of experimental MS peak intensities or areas (one column for each sample), with the first two columns representing \emph{m/z} and retention times of the peaks.

If not already available, the table can be obtained from an \code{xcmsSet} object (output of the \code{xcms} R package) through the \code{\link{table_xcms}} function.
}
  \item{info}{
Named list containing isotopic information, output of the \code{\link{isotopic_information}} function.
}
  \item{mass_shift}{
Maximum difference between theoretical and experimental mass. In other words, the mass accuracy.
}
  \item{RT}{
Expected retention time of the compound of interest.
}
  \item{RT_shift}{
Maximum difference between expected and experimental retention time of the peaks.
}
  \item{chrom_width}{
An estimate of the chromatographic peak width.
}
}


\value{
A matrix of extracted experimental isotopic patterns (one column for each sample), with the first two columns representing the exact \emph{m/z} and the retention times of the peaks.
}

\seealso{
\code{\link{table_xcms}} , \code{\link{isotopic_information}}
}