#' Compute theoretical isotopic patterns
#'
#' Function that computes the theoretical pattern of the specified compound, 
#' given the abundance of the labelling isotope X. 
#' The function returns a vector of normalised theoretical intensities, 
#' where the maximum signal is set to 100; 
#' the related masses are in the "target" vector contained in the input list.
#'
#' @param abundance Isotopic abundance of the labelling isotope X (either ^2H or ^13C); 
#' number between 0 and 1
#' @param info Named list containing isotopic information, 
#' output of the \code{\link{isotopic_information}} function
#' @param charge Natural number, denoting the charge state of the target adduct (1,2,3,...). If not provided, it is 1 by default 
#'
#' @return A vector representing the normalised isotopic pattern of the compound of interest, 
#' corresponding to the specified isotopic distribution
#' @export
#'
#' 
#' 
#' 
#' @author Ruggero Ferrazza
#' @references The function makes use of the \code{\link[ecipex]{ecipex}} R package
#' @seealso \code{\link{isotopic_information}}, \code{\link[ecipex]{ecipex}}
#' @keywords manip
#' 
#' 

pattern_from_abundance <-function(abundance, info, charge=1){
  
  # Get the table containing isotopic information
  isotopes <- info$isotopes
  
  # Modify the isotopes table with the abundance specified by the user
  isotopes[which(isotopes$element=="X")[1],"abundance"] <- 1 - abundance
  isotopes[which(isotopes$element=="X")[2],"abundance"] <- abundance
  
  theoret_pattern <- as.matrix(ecipex(info$compound, 
                                      isoinfo = isotopes, 
                                      limit = 1e-12, 
                                      id = FALSE, 
                                      sortby = "mass")[[1]])
  
  # Correct masses for charge state
  theoret_pattern[,1] <- theoret_pattern[,1]/charge
  
  
  # Group together signals coming from isotopologues with the same nucleon number,
  # and assign them the proper position
  theoretical_pattern <- unlist(lapply(info$target[-c(1,2)], function(x){  
    ind <- which(abs(x - theoret_pattern[,1]) <0.2/charge)
    return(sum(theoret_pattern[ind,2]))
                            }))
  # Normalise the theoretical pattern
  theoretical_pattern <- theoretical_pattern/max(theoretical_pattern)*100
  
  return(theoretical_pattern)
}
