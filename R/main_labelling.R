#' Main function of the package.
#' 
#' It computes the estimated X abundances of each sample, 
#' returning an object of the class "labelling".
#' 
#' @param peak_table A data.frame containing the integrated signals for the samples. 
#' @param compound the chemical formula of the compound of interest.
#' @param labelling character, either "H" or "C", specifying which is the labelling element.
#' @param mass_shift maximum shift allowed in the mass range.
#' @param RT expected retention time of the compund of interest.
#' @param RT_shift maximum shift allowed in the retention time range
#' @param chrom_width chromatographic width of the peaks
#' @param initial_abundance initial estimate for the abundance of the heaviest X
#' 
#' @details 
#' \itemize{
#'  \item{peak_table:}{ The first two columns of \code{peak_table} represent the mass and the retention time of the peaks; 
#'  the other columns represent peak intensities for each sample. 
#'  The table can be obtained using the function \code{table_xcms}}
#'  \item{compound:}{ Character vector, an X have to represent the elment with isotopic distribution to be fitted}
#'  \item{initial_abundance:}{ Numeric vector of the same length as the number of samples, with the initial estimate 
#'   for the abundance of the heaviest X isotope (either 2H or 13C). If provided, number between 0 and 100. Otherwise, NA}
#' }
#' 
#' @return An object of the class "labelling", which is a list containing the results from the fitting procedure:
#'  \item{compound}{Character vector specifying the chemical formula of the compound of interest, 
#'   with X being the element with unknown isotopic distribution (to be fitted)}
#'  \item{Best_estimate}{ Numeric vector representing the best estimated abundance 
#'  of the heaviest X isotope (either 2H or 13C). Number between 0 and 100.}
#'  \item{std_error}{ Numeric vector containing the standard errors of the estimates.}
#'  \item{dev_percent}{ The percentage deviations of the fitted theoretical patterns to the provided experimental patterns.}
#'  \item{x_scale}{ Vector containing the m/z signals of the isotopic patterns.}
#'  \item{y_exp}{ Matrix containing normalised experimental patterns, 
#'  where for each sample the most intense signal is set to 100.}
#'  \item{y_theor}{ Matrix of normalised fitted theoretical pattern (most intense signal set to 100 for each sample)}
#'  \item{warnings}{ Character vector containing possible warnings coming from the fitting procedure}
#'  
#'  
#' @examples
#' ## Examples are needed
#' @author Ruggero Ferrazza
#' @keywords manip



main_labelling <- function(peak_table, 
                           compound, 
                           labelling, 
                           mass_shift, 
                           RT, 
                           RT_shift, 
                           chrom_width, 
                           initial_abundance=NA){


  # Get some useful isotopic information, to be used in the coming functions
  info <- isotopic_information(compound, labelling)
  
  # Extract one pattern for each sample (each column of the peak_table data frame)
  experimental_patterns <- isotopic_pattern(peak_table, info, mass_shift, RT, RT_shift, chrom_width)
  
  # For each extracted pattern, find the X isotopic distribution that better fits the experimental data 
  fitted_abundances <- find_abundance(patterns=experimental_patterns, info, initial_abundance)
  
  return(fitted_abundances)
  
}
