\name{isotopic_pattern}
\alias{isotopic_pattern}

\title{
Extracting experimental isotopic patterns from a table of MS peaks
}

\description{
Function that extracts the experimental isotopic patterns of a specified compound from a table of peak signal intensities or areas.
}

\usage{
isotopic_pattern(peak_table, info, mass_shift, RT, RT_shift, chrom_width)
}


\arguments{
  \item{peak_table}{
Data frame containing experimental peak intensities or areas relative to a set of samples (one for each column), with the first two columns representing the m/z and the retention time of the peaks.

If not already available, the table can be obtained from an \code{xcmsSet} object (from the \code{xcms} R package) through the \code{\link{table_xcms}} function.
}
  \item{info}{
Named list containing isotopic information, output of the \code{\link{isotopic_information}} function.
}
  \item{mass_shift}{
Maximum mass shift between theoretical and experimental mass. In other words, the mass accuracy.
}
  \item{RT}{
Expected retention time of the compound of interest (in sec).
}
  \item{RT_shift}{
Maximum time shift between expected and experimental retention time of the peaks (in sec).
}
  \item{chrom_width}{
An estimate of the peak chromatographic width (in sec). 
}
}


\value{
A matrix of extracted experimental isotopic patterns: in the first two columns are the exact masses and the retention times of the peaks, and in the others are the extracted peak intensities or areas (one column for each sample).
}

\seealso{
\code{\link{table_xcms}} , \code{\link{isotopic_information}}
}