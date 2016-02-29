\name{isotopic_information}
\alias{isotopic_information}

\title{
Get useful isotopic information
}

\description{
This function gathers essential isotopic information required by the other functions of the \code{\link{IsotopicLabeling}} package.
}

\usage{
isotopic_information(compound, labeling)
}

\arguments{
  \item{compound}{
Character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted). 
}
  \item{labeling}{
Character, either "H" or "C", specifying the labeling element. 
}
}


\details{
The specified compound is not the neutral molecular species of interest, but the adduct observed by ESI-MS (such as protonated or sodiated species). 

In the chemical formula, the element with unknown abundance should be denoted by X. For example, the proton adduct of TAG 52:2, C55H103O6, should be written X55H103O6 for \eqn{^{13}}C labeling experiments, and C55X102HO6 for \eqn{^2}H labeling experiments. Note that in this last case only 102 hydrogen atoms have unknown isotopic distribution, since the one giving rise to the adduct comes from the solvent, and is considered to have fixed natural abundance.
}


\value{
A named list containing the following information:
  \item{compound}{The same as input.}
  \item{isotopes}{Table containing the natural isotopic abundances of the elements present in compound (numbers between 0 and 1). The two isotopes of element X are given NA values.}
  \item{target}{Named vector with the exact masses of all the possible isotopologues arising from the labeling isotope. M+0 is the monoisotopic mass (sum of the masses of the atoms using the lightest isotope for each element, X included); in M+1 one light isotope is replaced by its heaviest counterpart, and so forth.}
  \item{nX}{The number of X atoms. In other words, the number of atoms with unknown isotopic distribution.}
  \item{nTOT}{The total number of atoms of the labeling element (either H+X or C+X).}
}