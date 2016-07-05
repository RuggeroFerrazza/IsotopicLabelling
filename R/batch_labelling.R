#' Batch process a set of target analytes
#' 
#' This function batch processes LC- or GC-MS data, analyzing the isotopic patterns of a set of target analytes specified by the user. As with other functions of the package, it also requires chromatogrpahic information.
#' 
#' 
#' @param targets A data frame containing information on the target analytes. Columns 1 to 7 are: "name_compound" (the name of the target analytes), "compound", "labelling", "RT", "RT_shift", "chrom_width", "mass_shift" (the same parameters for \code{\link{main_labelling}} function). Columns from 8 onwards contain the initial estimates for the label abundances (one estimate for each sample). If not known, a single column of NA values can be entered
#' @param groups A factor containing the name of the group of each sample analysed; 
#' The function will calculate summary statistics for the samples belonging to the same group
#' @param plot_patterns,plot_residuals,plot_results Whether or not to plot the patterns, the residuals and a summary of the results. If so, pdf files are created in the working directory
#' @param save_results Whether to save the results of the estimates. If so, *.csv files are generated
#' 
#' @return batch_grouped_estimates A list as long as the number of target analytes, containing the summary of the fitted results (group estimates)
#' 
#' @export
#' 
#' @examples
#' # Get the sample dataset
#' data("xcms_obj")
#' 
#' # Convert the MS data set
#' peak_table <- table_xcms(xcms_obj)
#' 
#' # Get the example data frame containing target abalytes
#' data("targets")
#' 
#'  # Batch process the data
#'  batch_grouped_estimates <- batch_labelling(targets=targets, 
#'  groups=factor(c(rep("C12",4), rep("C13",4))),
#'  plot_patterns=FALSE, plot_residuals=FALSE, plot_results=FALSE, save_results=FALSE)
#'    
#' @author Ruggero Ferrazza
#' 
#' @seealso \link{main_labelling}, \link{group_labelling}, \link{save_labelling}, \link{plot.labelling}

batch_labelling <- function(targets, groups, plot_patterns=T, plot_residuals=F, plot_results=F, save_results=F){
  
  # attach targets data frame
  attach(targets)
  
  
  # Batch process
  batch_grouped_estimates <- list()
  
  for (i in 1:length(compound)){
    
    batch_fitted <- main_labelling(peak_table, compound=compound[i], labelling=labelling[i],
                                   mass_shift=mass_shift[i], RT=RT[i], RT_shift=RT_shift[i],
                                   chrom_width=chrom_width[i], initial_abundance=as.numeric(targets[i,8:ncol(targets)]))
    
    # plot the patterns
    if (plot_patterns) plot(x=batch_fitted, type="patterns", saveplots=T)
    
    # plot the residuals
    if (plot_residuals) plot(x=batch_fitted, type="residuals", saveplots=T)
    
    # plot the overall results
    if (plot_results) plot(x=batch_fitted, type="summary", saveplots=T)
    
    # save the results to a *.csv file
    if (save_results) save_labelling(batch_fitted)
    
    # Group the samples and obtain grouped estimates
    batch_grouped_estimates[[i]] <- group_labelling(batch_fitted, groups=groups)
    
  }
  
  names(batch_grouped_estimates) <- name_compound
  
  detach(targets)
  
  return(batch_grouped_estimates)
  
}