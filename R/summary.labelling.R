summary.labelling <-
function(fitted_abundances){
  # Function that produces a summary of the results from an object of class "labelling"
  
  # INPUT:
  # fitted_abundances: object of class "labelling"
  
  # OUTPUT:
  # results: matrix containing a summary of the fitted results. It has two rows, the first containing the estimated X abundances, and the second one containing the related standard errors coming from the fitting procedure 
  
  best_est <- fitted_abundances$best_estimate
  std_error <- fitted_abundances$std_error
  
  results <- rbind(best_est, std_error); row.names(results) <- c("Best Estimate [%]", "Standard Error [%]")
  
  return(results)
}
