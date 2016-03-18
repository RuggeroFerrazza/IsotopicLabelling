find_abundance <-
function(patterns, info, initial_abundance=NA){
  
  # Function that finds the isotopic distribution for element X that best reproduces the experimental patterns.
  
  # INPUT:
  # patterns: matrix containing the experimental patterns in its columns, with the first two representing the mass and the retention time of the peaks
  # info: named list, output of the "isotopic_information" function, containing useful information about the compound of interest and the isotopic distributions of its elements.
  # initial_abundance: numeric vector with length equal to the number of samples, with the initial estimate for the abundance of the heaviest X isotope (either 2H or 13C). If provided, number between 0 and 100. 
  
  # OUTPUT:
  # An object of the class "labelling", which is a list containing the results from the fitting procedure:
  # $compound: character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted)
  # $best_estimate: numeric vector representing the best estimated abundance of the heaviest X isotope (either 2H or 13C). Number between 0 and 100.
  # $std_error: numeric vector containing the standard errors of the estimates.
  # $dev_percent: the percentage deviations of the fitted theoretical patterns to the provided experimental patterns.
  # $x_scale: vector containing the m/z signals of the isotopic patterns.
  # $y_exp: matrix containing normalised experimental patterns, where for each sample the most intense signal is set to 100.
  # $y_theor: matrix of normalised fitted theoretical pattern (most intense signal set to 100 for each sample).
  # warnings: character vector containing possible warnings coming from the fitting procedure.
  
  #######      -------      #######

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
