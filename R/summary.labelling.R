#' Summary of the labelling fitting 
#'
#' Function that produces a summary of the results from an object of class \code{labelling}.
#'
#' @param fitted_abundances Object of class \code{labeling},
#' output of either \code{\link{main_labelling}} or \code{\link{find_abundance}} functions
#'
#' @return 
#' \item{results}{Matrix containing a summary of the fitted results. It has two rows, 
#' the first containing the estimated percentage isotopic abundances of the labelling isotope X (^2H or ^13C), 
#' and the second one containing the standard errors from the fitting procedure}
#' @export 
#'
#' @examples
#' 
#' @author Ruggero Ferrazza
#' 
#' 
#' @keywords manip
#' 
summary.labelling <- function(fitted_abundances){
  
  best_est <- fitted_abundances$best_estimate
  std_error <- fitted_abundances$std_error
  
  results <- rbind(best_est, std_error); 
  row.names(results) <- c("Best Estimate [%]", "Standard Error [%]")
  
  return(results)
}
