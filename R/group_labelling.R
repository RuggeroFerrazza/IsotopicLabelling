group_labelling <- 
function(fitted_abundances, groups){
  # Function that groups the labelling results in order to give a single estimated abundance for each sample group, with related standard error of the mean that takes into account both the standard error relative to each estimate coming from the fitting procedure, and the biological variability across samples
  
  # INPUT:
  # fitted_abundances: object of class "labelling"
  # groups: a factor containing the name of the group of each sample analysed; samples belonging to the same group will be put together.  
  
  # OUTPUT:
  # grouped_estimates: a data frame containing the results of the grouping. For each row (a group), it details: the number of samples belonging there, the averaged estimated percentage isotopic abundance of the labelling isotope, the standard error of the mean, the critical value for a 95% confidence interval of the t distribution with N-1 degrees of freedom, the lower and the upper 95% CI values.   
  
  ######  ------  ######
  
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