\name{isotopic_information}
\alias{isotopic_information}

\title{
Getting useful isotopic information
}

\description{
This function performs a pre-processing step to get some useful information about all the possible isotopic variants of the compound of interest. The gathered information will be used by other functions of the \code{\link{IsotopicLabeling}} package.
}

\usage{
isotopic_information(compound, labeling)
}

\arguments{
  \item{compound}{
Character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted). 
}
  \item{labeling}{
Character, either "H" or "C", specifying which is the labeling element. 
}
}


\details{
The compound to be specified should not be the chemical formula of the neutral molecular species of interest, but the adduct seen by ESI-MS (for example, protonated or sodiated species). 

As for the chemical formula, the element whose abundance is unknown should be denoted by X. Thus, if the interest is in the proton adduct of the TAG 52:2 species, C55H103O6, this should be written X55H103O6 for \eqn{^{13}}C labeling experiments, and C55X102HO6 for \eqn{^2}H labeling experiments. Note that in this last case only 102 hydrogen atoms are specified as having unkown isotopic distribution, whereas the one left out (coming from the solvents) will be considered to have fixed natural abundance.
}


\value{
A named list containing the following information:
  \item{compound}{The same as input.}
  \item{isotopes}{Table containing the natural isotopic abundances of the elements present in compound. The two isotopes of element X are given NA values.}
  \item{target}{Named vector containing the exact masses of all the possible isotopic variants of the species of interest. M+0 is the species where all the atoms (including X) are the lighest isotopes; in M+1 one light isotope is replaced by its heaviest counterpart, and so forth.}
  \item{nX}{The number of X atoms. In other words, the number of atoms with unkwown isotopic distribution.}
  \item{nTOT}{The total number of atoms of the labeling element (either H+X or C+X).}
}
