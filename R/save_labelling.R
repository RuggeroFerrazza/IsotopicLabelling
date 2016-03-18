save_labelling <-
function(fitted_abundances, path=getwd()){
  # Function that saves the obtained results to a csv file
  
  # INPUT:
  # fitted_abundances: object of class "labelling"
  # path: the directory where to save the csv file. If not specified, the results are saved in the working directory
  
  # OUTPUT:
  # The "COMPOUND_Estimated_Abundances.csv" file, containing the results
  
  #####  -----  #####
  table <- cbind(fitted_abundances$best_estimate, fitted_abundances$std_error, fitted_abundances$dev_percent, fitted_abundances$warnings)
  colnames(table) <- c("Best estimate [%]", "Standard Error [%]", "Percentage deviation [%]", "Fitting outcome messages/Warnings")
  
  filename <- paste(fitted_abundances$compound, "_Estimated_Abundances.csv", sep="")
  file <- paste(path, filename, sep="/")
  write.csv(x = table, file=file) 
}
