#' Compute a single estimated isotopic abundance for each sample group
#' 
#' This function groups the fitted abundances in order to give a single estimated value 
#' for each sample group, with related standard error of the mean that takes into account 
#' both the errors relative to each estimate from the fitting procedure, 
#' and the variability across samples.
#'
#' @param fitted_abundances Object of class \code{labelling}. 
#' It contains the results of the isotopic pattern analysis
#' @param groups A factor containing the name of the group of each sample analysed; 
#' The function will calculate summary statistics for the samples belonging to the same group
#'
#' @return A data frame containing the summary statistics calculated groupwise. 
#' For each row (a group), it details:
#' \item{N}{The number of samples in that group}
#' \item{Mean}{The averaged estimated percentage isotopic abundance of the labelling isotope}
#' \item{SE mean}{The standard error of the mean}
#' \item{t_crit}{The critical value for a 95\% confidence interval 
#' of the t distribution with N-1 degrees of freedom}
#' \item{Lower 95\% CI}{The lower 95\% confidence interval value}
#' \item{Upper 95\% CI}{The upper 95\% confidence interval value}
#' 
#' @details For each group, the average is simply computed by considering that the obtained individual values are representative of the population. 
#' 
#' As for the standard deviations, they are obtained using the law of total variance: the overall variance in each group is the sum of two distinct contributions, the first one related to the uncertainties associated in each sample estimate, and the second one arising from the spread of the estimates (biological variability).
#' 
#' @export
#'
#' @author Ruggero Ferrazza
#' @keywords manip
#' 
#' 
group_labelling <- function(fitted_abundances, groups){

  # Extract the estimated percentage abundances and the std errors of the fit
  estimates <- fitted_abundances$best_estimate
  std_err_fit <- fitted_abundances$std_error
  
  # Remove NA's from the data
  ind_NA <- which(is.na(estimates) | is.na(std_err_fit))
  
  if (length(ind_NA) !=0) {
    estimates <- estimates[-ind_NA]
    std_err_fit <- std_err_fit[-ind_NA]
    groups <- groups[-ind_NA]
  }
  
  # Compute the average for each group
  avg <- tapply(estimates, groups, mean)
  
  # Compute the variance due to biological variability
  var_biol <- tapply(estimates, groups, var)
  
  # Compute the variance due to single-sample errors
  var_within <- tapply(std_err_fit, groups, function(x){1/length(x)*sum(x^2)})
  
  # Compute the total variance
  var_TOT <- var_biol + var_within
  
  # Compute the std error of the MEAN
  N <- tapply(groups,groups,length)
  std_MEAN <- sqrt( var_TOT / N )
  
  # Compute the 95% Confidence intervals
  width <- tapply(groups, groups, function(x){ qt(.975, df=length(x)-1) })
  
  Lower <- avg - width*std_MEAN
  
  Upper <- avg + width*std_MEAN
  
  # Provide the output
  grouped_estimates <- data.frame(N, avg, std_MEAN, width, Lower, Upper);
  names(grouped_estimates) <- c("N", "Mean", "SE mean", "t_crit", "Lower 95% CI", "Upper 95% CI")
  
  return(grouped_estimates)
  
}