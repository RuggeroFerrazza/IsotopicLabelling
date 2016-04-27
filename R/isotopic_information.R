#' Get useful isotopic information
#' 
#' This function gathers essential isotopic information required by the other functions of the 
#' \code{\link{IsotopicLabelling}} package.
#' 
#' @param compound Character vector specifying the chemical formula of the compound of interest, 
#' with X being the element with unknown isotopic distribution (to be fitted)
#' @param labelling Character, either "H" or "C", specifying the labelling element 
#' 
#' @return A list with the following elements:
#' \item{compound}{The same as input}
#' \item{target}{Named vector with the exact masses of all the possible isotopologues 
#' arising from the labelling isotope. 
#' M+0 is the monoisotopic mass (sum of the masses of the atoms using the lightest isotope for each element, X included); 
#' in M+1 one light isotope is replaced by its heaviest counterpart, and so forth}
#' \item{isotopes}{Table containing the natural isotopic abundances of the elements present in compound (numbers between 0 and 1).
#'  The two isotopes of element X are given NA value}
#' \item{nX}{The number of X atoms. In other words, the number of atoms with unknown isotopic distribution}
#' \item{nTOT}{The total number of atoms of the labelling element (either H+X or C+X)} 
#' 
#' @details  The specified compound is not the neutral molecular species of interest, 
#' but the adduct observed by ESI-MS (such as protonated or sodiated species). 
#' In the chemical formula, the element with unknown abundance should be denoted by X. 
#' For example, the proton adduct of TAG 52:2, C55H103O6, should be written X55H103O6 for 
#' ^13C labelling experiments, and C55X102HO6 for ^2H labelling experiments. 
#' Note that in this last case only 102 hydrogen atoms have unknown isotopic distribution, 
#' since the one giving rise to the adduct comes from the solvent, 
#' and is considered to have fixed natural abundance.
#' 
#' @export
#'
#' @examples
#' info <- isotopic_information(compound="X40H77NO8P", labelling="C") 
#' # This is the case for [PC 32:2+H]+ in a ^13C-labelling experiment 
#' @author Ruggero Ferrazza
#' @keywords manip
#' 


isotopic_information <- function(compound, labelling){

  # Check that labelling is correct
  if (labelling !="H" & labelling !="C") stop("Check the labelling character: it should be either H or C")
  
  
  # Introduce X in the isotopes data 
  isotopes <- nistiso[,1:3]
  X_new <- isotopes[which(isotopes$element==labelling),] 
  X_new$element <- "X"
  
  isotopes <- rbind(isotopes, X_new)
  
  
  # Compute the mass difference between heavier and lighter isotope
  mass_diff <- abs(diff(X_new[,"mass"]))
  
  
  # In isotopes data frame keep only the elements of interest
  DF <- strapply(compound, 
                 "([A-Z][a-z]*)(\\d*)", 
                 ~ c(..1, if (nchar(..2)) ..2 else 1), 
                 simplify = ~ as.data.frame(t(matrix(..1, 2)), stringsAsFactors = FALSE)) 
  DF[[2]] <- as.numeric(DF[[2]]) 
  

  isotopes <- isotopes[isotopes$element %in% DF[,1],] 
  
  
  # Compute the exact mass of the species of interest (suppose that X has natural abundance)
  DF[[3]] <- apply(DF, 1, function(x){a <- which(x[1]==isotopes$element)
                                      a <- a[which.max(isotopes$abundance[a])]
                                      return(isotopes$mass[a])})
  
  exact_mass <- sum(DF[[2]]*DF[[3]])
  
  
  # Set the X abundance to be unknown
  isotopes$abundance[which(isotopes$element == "X")] <- NA
  row.names(isotopes) <- seq(nrow(isotopes))
  
  # Extract total number of atoms of the element being labelled
  nTOT <- DF[which(DF[,1]==labelling),2]
  if (length(nTOT)==0) nTOT <- 0
  
  nX <- DF[which(DF[,1]=="X"),2]
  if (length(nX)==0) nX <- 0
  
  nTOT <- nTOT + nX
  
  
  # Create the target vector, containing the exact masses of all the possible isotopic variants arising from X
  # Lowest mass: 2 mass units below the monoisotopic mass
  # Highest mass: 2 mass units above the mass corresponding to all X atoms having been labelled with the heaviest isotope
  
  target <- round(seq(from=exact_mass-2*mass_diff, by=mass_diff, length=nX+5), digits=4)
  names(target) <- c("M-2", "M-1", paste("M+", 0:(nX+2), sep=""))
  
  return(list(compound=compound, isotopes=isotopes, target=target, nX=nX, nTOT=nTOT))
  
}
