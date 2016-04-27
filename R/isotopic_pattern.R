#' Extract experimental isotopic patterns from a table of MS peaks
#'
#' Function that extracts the experimental isotopic patterns of a specified compound from a data frame 
#' containing MS peak intensities or areas.
#'
#' @param peak_table Data frame of experimental MS peak intensities or areas (one column for each sample), 
#' with the first two columns representing \emph{m/z} and retention times of the peaks
#' @param info Named list containing isotopic information, 
#' output of the \code{\link{isotopic_information}} function
#' @param mass_shift Maximum difference between theoretical and experimental mass. 
#' In other words, the expected mass accuracy
#' @param RT Expected retention time of the compound of interest
#' @param RT_shift Maximum difference between expected and experimental retention time of the peaks
#' @param chrom_width An estimate of the chromatographic peak width
#'
#' @return  A matrix of extracted experimental isotopic patterns (one column for each sample), 
#' with the first two columns representing the exact \emph{m/z} and the retention times of the peaks
#' 
#' @details The table can be obtained from an \code{xcmsSet} 
#' object (output of the \code{xcms} R package) through the \code{\link{table_xcms}} function.
#' 
#' 
#' @export
#'
#' @examples
#' experimental_patterns <- isotopic_pattern(peak_table, info, mass_shift=0.05, 
#' RT=285, RT_shift=20, chrom_width=7)
#' 
#' @author Ruggero Ferrazza
#' @seealso \code{\link{table_xcms}} , \code{\link{isotopic_information}}
#' @keywords manip

isotopic_pattern <-function(peak_table, info, mass_shift, RT, RT_shift, chrom_width){
  
  tmp_list <- lapply(info$target, function(x){ind <- which( (abs(peak_table$mz - x) < mass_shift) & (peak_table$rt < (RT + RT_shift) ) & (peak_table$rt > (RT - RT_shift) ) )
                                  return(data.frame(ind=ind, rt=peak_table[ind,"rt"]))
                                 })
  
    # Extract the retention times of all the peaks
  rt_overall <- sort(unique(unlist(lapply(tmp_list, function(x){x$rt}), use.names=F)))
  
  rt_grouped <- apply(abs(outer(rt_overall,rt_overall,'-')), 2, function(u) list(rt_overall[u<=chrom_width]))
  rt_grouped <- unique(lapply(rt_grouped, "[[", 1))
  
  rt_candidates <- sapply(rt_grouped, mean)
  rt_best <- rt_candidates[which.min(abs(rt_candidates - RT))]

  # Define the matrix where to put the signals
  patterns <- matrix(0, nrow=length(info$target), ncol=(ncol(peak_table))) 
  row.names(patterns) <- names(info$target)
  colnames(patterns) <- colnames(peak_table)

  
  for (i in 1:length(info$target)){
    
    ind <- which( (abs(peak_table$mz - info$target[i]) < mass_shift) & (abs(peak_table$rt - rt_best)<= chrom_width) )
    
    if (length(ind)>=1){
      ind <- ind[which.min(abs(peak_table$rt[ind] - rt_best))]

      patterns[i,] <- as.numeric(peak_table[ind,])

    }
    
  }

  # Check that the most intense signals do not come from M-2 or M-1 (which could arise from the same species as the target, with one more unsaturation)
  # If so, "delete" the whole experimental pattern
  
  if (ncol(patterns)>=3){
    max_pos <- apply(patterns[,-c(1,2)], 2, which.max)
    
    patterns[,c(F,F,max_pos <=2)]  <- 0
  }
  
  patterns[is.na(patterns)] <- 0
  
  # Cut out the first two masses (M-2 and M-1)
  patterns <- patterns[-c(1,2),]

  # Check that each pattern has at least two signals different from 0, otherwise set all to 0
  ind_pat <- apply(patterns[,-c(1,2)], 2, function(c)sum(c!=0)) > 1
  patterns[,c(F,F,!ind_pat)] <- 0
  
  # Add the exact masses to the patterns matrix
  patterns[,"mz"] <- info$target[-c(1,2)]
  patterns[ which(patterns[,"rt"]==0), "rt"] <- NA
  
  # Return the obtained patterns
  return(patterns)
  
}
