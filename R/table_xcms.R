table_xcms <-
function(xcms_obj){
  # Function that properly convert an xcmsSet object, from package xcms, into a table of peaks
  
  # INPUT:
  # xcms_obj: an xcmsSet  object, coming from the processing of MS data with the xcms R package

  # OUTPUT:
  # peak_table: data frame containing the signals in the xcmsSet object at input, but with the proper format: it is a data frame with the first two columns representing mass and retention time of the related peaks
  
  
  ######  ------  ######

  # Check that the file in input is an xcmsSet object
  
  if (class(xcms_obj) != "xcmsSet") stop("ERROR: The provided object is not an xcmsSet object")
  
  peak_table <- peakTable(xcms_obj)
  
  n_files <- length(sampnames(xcms_obj))
  
  peak_table <- peak_table[,c(which(colnames(peak_table)=="mz" | colnames(peak_table)=="rt"),
                (ncol(peak_table)-n_files+1):ncol(peak_table))]
  
  return(peak_table)
}
