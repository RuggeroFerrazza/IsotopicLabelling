\name{isotopic_pattern}
\alias{isotopic_pattern}

\title{
Extracting experimental isotopic patterns from a table of MS peaks
}

\description{
Function that extracts the experimental isotopic patterns of a specificed compound from a table of peak signal intensities or areas.
}

\usage{
isotopic_pattern(peak_table, info, mass_shift, RT, RT_shift, chrom_width)
}


\arguments{
  \item{peak_table}{
Data frame containing experimental peak intensities or areas relative to a set of samples (one sample for each column), with the first two columns representing the m/z and the retention time of the peaks.

If not available, the table can be obtained from the function \code{\link{table_xcms}}, starting from an xcmsSet object (from the \code{xcms} R package).
}
  \item{info}{
Named list containing isotopic information, output of the \code{\link{isotopic_information}} function.
}
  \item{mass_shift}{
Maximum shift allowed in the mass range between theoretical and experimental mass. In other words, the accuracy of the measurements.
}
  \item{RT}{
Expected retention time of the compound of interest (in sec).
}
  \item{RT_shift}{
Maximum shift allowed for the peaks in the retention time range (in sec).
}
  \item{chrom_width}{
An estimate of the chromatographic width of the peaks (in sec). 
}
}


\value{
A matrix of extracted experimental isotopic patterns: in the first two columns are the exact masses and the retention times of the peaks, and in the others are the extracted peak intensities or areas (one column for each sample).
}

\seealso{
\code{\link{table_xcms}} , \code{\link{isotopic_information}}
}