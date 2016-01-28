main_labeling <-
function(peak_table, compound, labeling, mass_shift, RT, RT_shift, chrom_width, initial_abundance=NA){
  # Main function of the package. It computes the estimated X abundances of each sample, returning an object of the class "labeling".
  
  # INPUT:
  # peak_table: data frame containing the integrated signals for the samples. The first two columns represent the mass and the retention time of the peaks; the other columns represent peak intensities for each sample. The table can be obtained using the function "table_xcms"
  # compound: character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted)
  # labeling: character, either "H" or "C", specifying which is the labeling element
  # mass_shift: maximum shift allowed in the mass range
  # RT: expected retention time of the compund of interest
  # RT_shift: maximum shift allowed in the retention time range
  # chrom_width: chromatographic width of the peaks
  # initial_abundance: numeric vector of the same length as the number of samples, with the initial estimate for the abundance of the heaviest X isotope (either 2H or 13C). If provided, number between 0 and 1. Otherwise, NA

  # OUTPUT:
  # An object of the class "labeling", which is a list containing the results from the fitting procedure:
  # $compound: character vector specifying the chemical formula of the compound of interest, with X being the element with unknown isotopic distribution (to be fitted)
  # $Best_estimate: numeric vector representing the best estimated abundance of the heaviest X isotope (either 2H or 13C). Number between 0 and 1.
  # $std_error: numeric vector containing the standard errors of the estimates.
  # $dev_percent: the percentage deviations of the fitted theoretical patterns to the provided experimental patterns.
  # $x_scale: vector containing the m/z signals of the isotopic patterns.
  # $y_exp: matrix containing normalized experimental patterns, where for each sample the most intense signal is set to 100.
  # $y_theor: matrix of normalized fitted theoretical pattern (most intense signal set to 100 for each sample).
  # warnings: character vector containing possible warnings coming from the fitting procedure.
  
  
  #####  -----  #####
  

  # Get some useful isotopic information, to be used in the coming functions
  info <- isotopic_information(compound, labeling)
  
  # Extract one pattern for each sample (each column of the peak_table data frame)
  experimental_patterns <- isotopic_pattern(peak_table, info, mass_shift, RT, RT_shift, chrom_width)
  
  # For each extracted pattern, find the X isotopic distribution that better fits the experimental data 
  fitted_abundances <- find_abundance(patterns=experimental_patterns, info, initial_abundance)
  
  return(fitted_abundances)
  
}
