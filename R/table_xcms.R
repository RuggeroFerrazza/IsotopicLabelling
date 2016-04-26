#' Process \code{xcmsSet}  
#'
#' Function that properly convert an xcmsSet object, 
#' from package xcms, into a table of peaks
#'
#' @param xcms_obj an xcmsSet object
#' 
#' @return 
#' \item{peak_table}{data frame extracted from the \code{xcmsSet} object. 
#' The first two columns representing mass and retention 
#' time of the related peaks}
#'
#' @note The output data frame, required by other functions of the 
#' \code{\link{IsotopicLabelling}} R package, 
#' can be obtained in a number of other independent ways, 
#' such as through proprietary software of the vendor of the MS instrument. 
#'
#' @author Ruggero Ferrazza
#' 
#' @examples
#' data(xcms_obj)
#' peak_table <- table_xcms(xcms_obj) 
#' 
#' @keywords manip
#' @export

table_xcms <- function(xcms_obj){
  
  # Check that the file in input is an xcmsSet object
  
  if (class(xcms_obj) != "xcmsSet") stop("ERROR: The provided object is not an xcmsSet object")
  
  peak_table <- peakTable(xcms_obj)
  
  n_files <- length(sampnames(xcms_obj))
  
  peak_table <- peak_table[,c(which(colnames(peak_table)=="mz" | colnames(peak_table)=="rt"),
                (ncol(peak_table)-n_files+1):ncol(peak_table))]
  
  return(peak_table)
}
