plot.labeling <-
function(fitted_abundances, path=NA, what="patterns"){
  # Function that plots, for each sample, both its normalized experimental pattern and the best fitted theoretical pattern. The best estimated X abundance is also given 
  # In addition, another plot is produced containing the fitted abundances with related standard errors 
  # The plots are saved to the specified directory. If not defined, a folder is created in the working directory, with the name of the compound 
  
  # INPUT:
  # fitted_abundances: object of class "labeling"
  # path: the directory where to save the plots. If not specified, a folder with the name of the compound is created in the working directory, and the plots are saved there
  # what: character vector specifying what to plot. Either "patterns" (the default value) or "residuals". In this last case, the residuals are plotted both to screen and to file.
  
  # OUTPUT:
  # The resulting plots
  
  #####  -----  #####
  if (what=="patterns"){
  
      if (is.na(path)){ path <- paste(getwd(), fitted_abundances$compound, sep="/")
                    dir.create(path=path) }
  
  
      sample_name <- names(fitted_abundances$best_estimate)
  
      # Plot the results 

      for (k in 1:length(sample_name)){
    
        plot(fitted_abundances$x_scale, fitted_abundances$y_exp[,k], type="h", main=paste(sample_name[k], " ,   Compound: ",  fitted_abundances$compound), xlab="Target mass", ylab="Normalized intensity", ylim=c(0,110), cex.main=1)
        text=100*fitted_abundances$best_estimate[k]
    
        mtext(paste("Fitted X abundance: ", sprintf("%1.3f%%", text)))
    
        points(fitted_abundances$x_scale, fitted_abundances$y_theor[,k], col=2, pch=16, cex=.5)
        points(fitted_abundances$x_scale[1], -2.5, pch=17, col="blue")
    
        legend("top", legend=c("Experimental pattern", "Theoretical pattern"), lty=c(1,0), pch=c(0,16), pt.cex=c(0,0.6), col=c(1,2), cex=0.7, horiz=TRUE)
    
        dev.copy(png,filename=paste(paste(path, sample_name[k], sep="/"), ".png", sep=""), width=10, height=8, units="in", res=500)
    
        dev.off()
    
      }
  
      par(mar=c(6.4,4.1,4.1,2.1))
      plot(fitted_abundances$best_estimate, type="p", pch=16, cex=0.6, xaxt="n", main=paste("Estimated abundances", " ,   Compound: ",  fitted_abundances$compound), ylab="Fitted X abundance", xlab="", cex.main=1)
  
      axis(side=1, at=1:length(sample_name), labels=sample_name, las=2, cex.axis=0.7)
  
      segments(x0=1:length(sample_name), y0=fitted_abundances$best_estimate - fitted_abundances$std_error, x1=1:length(sample_name), y1=fitted_abundances$best_estimate + fitted_abundances$std_error, col="gray")
  
      dev.copy(png,filename=paste(paste(path, "# Estimated abundances", sep="/"), ".png", sep=""), width=10, height=8, units="in", res=500)
  
      dev.off()
      
  } else if (what=="residuals") {
    
    if (is.na(path)){ path <- paste(getwd(), fitted_abundances$compound, sep="/")
    dir.create(path=path) }
    
    
    sample_name <- names(fitted_abundances$best_estimate)
    
    # Plot the residuals
    
    for (k in 1:length(sample_name)){
      
      plot(fitted_abundances$x_scale, fitted_abundances$residuals[,k], type="h", main=paste("Residuals for ", sample_name[k]), xlab="Target mass", ylab="Normalized intensity", cex.main=1)
      text=100*fitted_abundances$best_estimate[k]
      
      mtext(paste("Fitted X abundance: ", sprintf("%1.3f%%", text)))
      
      abline(h=0, col="gray")
      
      dev.copy(png,filename=paste(paste(path, paste("Residuals_", sample_name[k], sep=""), sep="/"), ".png", sep=""), width=10, height=8, units="in", res=500)
      
      dev.off()
      
    }
    
  }
  
  
}
