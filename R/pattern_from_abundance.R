pattern_from_abundance <-
function(abundance, info){
  # Function that, given in input the isotopic abundance of either 2H or 13C, computes the theoretical pattern
  # The function returns a vector containing the normalised theoretical intensities (maximum signal set to 100) corresponding to the masses in the "target" vector 
  
  # INPUT:
  # abundance: isotopic abundance of either 2H or 13C
  # info: list containing useful information to be used. Output of the "isotopic_information" function
  
  # OUTPUT:
  # theoretical_pattern: vector containing the normalised intensities of the theoretical pattern expected if X had the isotopic abundance given at the input
  
  
  ###################

  # Get the table containing isotopic information
  isotopes <- info$isotopes
  
  # Modify the isotopes table with the abundance specified by the user
  isotopes[which(isotopes$element=="X")[1],"abundance"] <- 1 - abundance
  isotopes[which(isotopes$element=="X")[2],"abundance"] <- abundance
  
  theoret_pattern <- as.matrix(ecipex(info$compound, isoinfo = isotopes, limit = 1e-12, id = FALSE, sortby = "mass")[[1]])
  
  # Group together signals coming from isotopologues with the same nucleon number, and assign them the proper position
  theoretical_pattern <- unlist(lapply(info$target[-c(1,2)], function(x){  ind <- which(abs(x - theoret_pattern[,1]) <0.2)
                               return(sum(theoret_pattern[ind,2]))
                            }))
  
 
  # Normalise the theoretical pattern
  theoretical_pattern <- theoretical_pattern/max(theoretical_pattern)*100
  
  return(theoretical_pattern)
  
}
