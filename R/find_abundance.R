#' Fit experimental isotopic patterns
#' 
#' Function that takes each of the provided experimental MS isotopic patterns, 
#' and fits the best theoretical pattern that reproduces it through a weighted non-linear least squares 
#' procedure.
#'
#'
#' @param patterns A matrix of experimental isotopic patterns (one column for each sample), 
#' with the first two columns representing \emph{m/z} and retention time of the corresponding peaks
#' @param info Named list containing isotopic information, output of the \code{\link{isotopic_information}} function
#' @param initial_abundance Either NA, or a numeric vector of length equal to the number of samples, 
#' with the initial guesses on the percentage isotopic abundance of the labelling isotope 
#' (denoted as X, it can be either ^2H or ^13C). If provided, numbers between 0 and 100
#'
#' @return An object of class \code{labelling}, 
#' which is a list containing the results of the fitting procedure:
#' \item{compound}{Character vector specifying the chemical formula of the compound of interest, 
#' with X being the element with unknown isotopic distribution (to be fitted)}
#' \item{best_estimate}{Numeric vector of length equal to the number of samples, 
#' containing the estimated percentage abundances of the labelling isotope X 
#' (either ^2H or ^13C). Numbers between 0 and 100}
#' \item{std_error}{Numeric vector with the standard errors of the estimates, 
#' output of the \code{nls} fitting procedure}
#' \item{dev_percent}{Numeric vector with the percentage deviations between best fitted and related experimental patterns}
#' \item{x_scale}{Numeric vector containing the \emph{m/z} values relative to the signals of the experimental patterns}
#' \item{y_exp}{Matrix of normalised experimental isotopic patterns (one column for each sample). 
#' The most intense signal of each pattern is set to 100}
#' \item{y_theor}{Matrix of normalised fitted theoretical isotopic patterns (one column for each sample). 
#' The most intense signal of each pattern is set to 100}
#' \item{residuals}{Matrix of residuals: each column is the difference between experimental and best fitted theoretical patterns}
#' \item{warnings}{Character vector with possible warnings from the \code{nls} fitting procedure}
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' fitted_abundances <- find_abundance(patterns, info, initial_abundance=NA)
#' }
#' 
#' @author Ruggero Ferrazza
#' @seealso \code{\link{isotopic_information}}

find_abundance <- function(patterns, info, initial_abundance=NA){
  tmp_results <- list()
  
  
      analysis_X <- function(pattern, info, initial_ab=NA){
            
            # Create a vector of masses
            target <- info$target[-c(1,2)]
            
            # Create the list to return in case of errors
            error_list <- list(compound=info$compound, best_estimate=NA, std_error=NA, dev_percent=NA, x_scale=target, y_exp=pattern/max(pattern)*100, y_theor=rep(NA, times=length(pattern)), residuals=rep(NA, times=length(pattern)), warnings="An error occurred")
    
            if (sum(pattern)==0) return(error_list)
    
            # Normalise the experimental pattern (max intensity = 100)
            pattern <- pattern/max(pattern)*100
    
            # Find and store the mass of the most intense signal
            mass_max <- target[which.max(pattern)]
    
            # First, rough estimate of the X abundance (either 2H or 13C), using mass_max and the exact mass 
            # If the user inserts a first estimate (initial_abundance), skip this step
            if (is.na(initial_ab)) initial_ab <-  unname(round((mass_max - target[1])/info$nX, digits=3)) 
            if (initial_ab <0) initial_ab <- 0 
            if (initial_ab >1) initial_ab <- 1 
    
            # Fitting procedure to find the best estimate for the X isotopic abundance 
            # The signals of the pattern are given weights proportional to the square root of their intensity, so as to give less importance to noise
            # Define a function of only one parameter, the abundance, which will have to be fitted by the nls function
            pattern_fit <- function(abundance) { pattern_from_abundance(abundance, info=info)}
    
            fit <- nls(formula=pattern~pattern_fit(abundance), start=list(abundance=initial_ab), control=list(maxiter=50, tol=5e-8, warnOnly=T), algorithm="port", weights=sqrt(pattern), na.action=na.exclude, lower=0, upper=1)
    
            if (inherits(try(summary(fit), silent=TRUE),"try-error")) return(error_list)
    
            warnings <- fit$convInfo$stopMessage
    
            return(list(best_estimate=summary(fit)$coefficients[1]*100, std_error=summary(fit)$coefficients[2]*100, dev_percent=(sqrt(sum((summary(fit)$residuals)^2)/ sum(pattern^2))*100), x_scale=target, y_exp=pattern, y_theor=pattern_fit(summary(fit)$coefficients[1]), residuals=pattern-pattern_fit(summary(fit)$coefficients[1]), warnings=warnings))
    
      }
  
  
  for (i in 3:ncol(patterns)){
    
    tmp_results[[i-2]] <- analysis_X(pattern=patterns[,i], info=info, initial_ab=initial_abundance[i-2]/100)

  }
  names(tmp_results) <- colnames(patterns[,-c(1,2)])
      
  # Create the output list from the results obtained 
  best_estimate <- sapply(tmp_results, "[[", "best_estimate")     
  std_error <- sapply(tmp_results, "[[", "std_error")  
  dev_percent <- sapply(tmp_results, "[[", "dev_percent")
  x_scale <- patterns[,"mz"]
  y_exp <- sapply(tmp_results, "[[", "y_exp")
  y_theor <- sapply(tmp_results, "[[", "y_theor")
  residuals <- sapply(tmp_results, "[[", "residuals")
  warnings <- sapply(tmp_results, "[[", "warnings")

  results <- list(compound=info$compound, best_estimate=best_estimate, std_error=std_error, dev_percent=dev_percent, x_scale=x_scale, y_exp=y_exp, y_theor=y_theor, residuals=residuals, warnings=warnings)
  class(results) <- "labelling"
  
  return(results)

}
