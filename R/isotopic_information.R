isotopic_information <-
function(compound, labelling){
  
  # Function that provides a list of isotopic information on the compound at input
  
  # INPUT:
  # compound: character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted)
  # labelling: character, either "H" or "C", specifying which is the labelling element
  
  # OUTPUT:
  # a named list containing the following information:
  # compound: the same as input
  # isotopes: table containing the natural isotopic abundances of the elements present in compound. The X element is given NA values
  # target: vector containing the exact masses of all the possible isotopic variants of the species of interest
  # nX : the number of X atoms
  # nTOT: the total number of atoms of the labelled element (either H+X or C+X)
  
  ##### 

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
