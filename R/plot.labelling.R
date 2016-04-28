#' Plot method for \code{labelling} objects
#' 
#' Produces different types of summary plots for a \code{labelling} object.
#'
#' @param x Object of class \code{labelling}
#' @param type The type of output produced. Available options are "patterns", "residuals", "summary"
#' @param saveplots Should the plots be saved as a PDF file?
#'
#' @return One or more plots
#' @details The default (type 'patterns') plot shows, for each sample in the class \code{labelling} object, 
#' the normalized experimental pattern superimposed to its fitted theoretical pattern.
#' By setting type to 'residuals', the function plots the residuals (the differences between experimental and best fitted theoretical patterns).
#' Type 'summary' produces a summary plot showing the estimated percentage abundances with related standard errors.
#' If \code{saveplot} is TRUE, the plots are saved to a PDF file in the working directory.
#' @export
#'
#' @examples
#' ## to be added
#' plot(x=fitted_abundances, type="patterns", saveplots=TRUE) 
#' plot(x=fitted_abundances, type="residuals", saveplots=TRUE) 
#' plot(x=fitted_abundances, type="summary", saveplots=TRUE)
#' 
#' 
#' @author Ruggero Ferrazza
#' @keywords hplot

  

plot.labelling <- function(x, type="patterns", saveplots=F, ...){
  fitted_abundances <- x
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
        text=paste("Fitted X Abundance: (", sprintf("%1.3f", fitted_abundances$best_estimate[k]), "+/-", sprintf("%1.3f", fitted_abundances$std_error[k]), ")\ %")
    
        mtext(text, cex=0.8)
    
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
      text=paste("Fitted X Abundance: (", sprintf("%1.3f", fitted_abundances$best_estimate[k]), "+/-", sprintf("%1.3f", fitted_abundances$std_error[k]), ")\ %")
      
      mtext(text, cex=0.8)
      
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
