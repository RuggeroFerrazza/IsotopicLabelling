#' Export to csv
#'
#' Function that saves the obtained results to a csv file.
#'
#' @param fitted_abundances Object of class \code{labelling}
#' @param path The directory where to save the csv file. 
#' If not specified, the results are saved in the working directory
#'
#' @return The "COMPOUND_Estimated_Abundances.csv" file, containing the results of the analysis. 
#' For each sample (one for each row) there are four columns:
#' \enumerate{
#' \item The estimated percentage abundance of the labelling isotope (either ^2H or ^13C);
#' \item The related standard error coming from the \code{nls} fitting procedure;
#' \item The percentage deviation between theoretical and experimental isotopic patterns;
#' \item The outcome message from the fitting procedure, to undersand whether there have been any convergence problems.
#' }
#' @export
#' @author Ruggero Ferrazza
#' @examples
#' 
#' @keywords IO
#' @seealso \link{main_labelling}
#' 
save_labelling <-function(fitted_abundances, path=getwd()){
  table <- cbind(fitted_abundances$best_estimate, 
                 fitted_abundances$std_error, 
                 fitted_abundances$dev_percent,
                 fitted_abundances$warnings)
  colnames(table) <- c("Best estimate [%]", "Standard Error [%]", 
                       "Percentage deviation [%]", "Fitting outcome messages/Warnings")
  
  filename <- paste(fitted_abundances$compound, "_Estimated_Abundances.csv", sep="")
  file <- paste(path, filename, sep="/")
  write.csv(x = table, file=file) 
}

