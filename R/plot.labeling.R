plot.labeling <- function(fitted_abundances, type="patterns", saveplots=F){
  # Function to plots the results of the Isotopic Labeling analysis
  # Three types of plots can be produced, depending on the parameters "type":
  # 1. The experimental patterns together with the best fitted theoretical patterns (default type of plot)
  # 2. The residuals 
  # 3. A summary of the estimated abundances
  
  # INPUT:
  # fitted_abundances: object of class "labeling"
  # type: character vector specifying what to plot: "patterns" (the default value), "residuals" or "summary"
  # saveplots: logical. If TRUE, the plots are saved to a *.pdf file in the working directory.
   
  # OUTPUT:
  # The resulting plots
  
  #####  -----  #####
  
  # Save default parameters
  plot.new()
  old.par <- par(no.readonly = T)
  
  sample_name <- names(fitted_abundances$best_estimate)
  
    if (type=="patterns"){

      # Plot the results 

      if (saveplots==T) { 
        pdf(paste(fitted_abundances$compound, "_Isotopic_Patterns", ".pdf", sep=""), width=6.5, height=3) 
        par(mar=c(3,3,2.5,0.1), mgp=c(2,0.6,0))
        } else par(mar=c(4,4.5,4,0.1), mfrow=c(3,1), ask=T)
      
      for (k in 1:length(sample_name)){
    
        plot(fitted_abundances$x_scale, fitted_abundances$y_exp[,k], type="h", main=paste(sample_name[k], " ,   Compound: ",  fitted_abundances$compound), xlab="Target mass", ylab="Normalised intensity", ylim=c(0,110), cex.main=1)
        text=fitted_abundances$best_estimate[k]
    
        mtext(paste("Fitted X abundance: ", sprintf("%1.3f%%", text)), cex=0.8)
    
        points(fitted_abundances$x_scale, fitted_abundances$y_theor[,k], col=2, pch=16, cex=.5)
        points(fitted_abundances$x_scale[1], -2.5, pch=17, col="blue")
    
        legend("top", legend=c("Experimental pattern", "Theoretical pattern"), lty=c(1,0), pch=c(0,16), pt.cex=c(0,0.6), col=c(1,2), cex=0.7, horiz=TRUE)
    
      }
  
     
  } else if (type=="residuals") {
    

    # Plot the residuals
    
    if (saveplots==T) {
      pdf(paste(fitted_abundances$compound, "_Residuals", ".pdf", sep=""), width=6.5, height=3)
      par(mar=c(3,3,2.5,0.1), mgp=c(2,0.6,0))
      } else par(mar=c(4,4.5,4,0.1), mfrow=c(3,1), ask=T)


    for (k in 1:length(sample_name)){
      
      plot(fitted_abundances$x_scale, fitted_abundances$residuals[,k], type="h", main=paste("Residuals for ", sample_name[k]), xlab="Target mass", ylab="Normalised intensity", cex.main=1)
      text=fitted_abundances$best_estimate[k]
      
      mtext(paste("Fitted X abundance: ", sprintf("%1.3f%%", text)), cex=0.8)
      
      abline(h=0, col="gray")
      
    }
    
  } else if (type=="summary"){
    
    if (saveplots==T) pdf(paste(fitted_abundances$compound, "_Summary", ".pdf", sep=""), width=10, height=7) 
    par(mar=c(6.4,4.1,4.1,2.1))
    
    plot(fitted_abundances$best_estimate, type="p", pch=16, cex=0.6, xaxt="n", main="Summary of the Estimated Abundances", ylab="Fitted X abundance  [%]", xlab="", cex.main=1.2)
    mtext(paste("Compound: ", fitted_abundances$compound))
    axis(side=1, at=1:length(sample_name), labels=sample_name, las=2, cex.axis=0.7)
    
    segments(x0=1:length(sample_name), y0=fitted_abundances$best_estimate - fitted_abundances$std_error, x1=1:length(sample_name), y1=fitted_abundances$best_estimate + fitted_abundances$std_error, col="gray")
    
  }
  
  if (saveplots==T) dev.off() else par(old.par)
}
