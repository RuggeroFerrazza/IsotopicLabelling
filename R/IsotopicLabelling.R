#' IsotopicLabelling-package
#'
#' The \code{IsotopicLabelling} package allows to analyse the isotopic patterns in MS data 
#' obtained in isotopic labelling experiments.  From the experimental patterns, 
#' the package estimates the isotopic abundance of the stable isotope employed in 
#' the labelling experiment (either ^2H or ^13C) inside a specified compound. 
#' 
#' 
#' @section Details:
#'Given a data frame of LC-MS or GC-MS peak intensities or areas (one column for each sample to analyse), 
#' the \code{\link{IsotopicLabelling}} package first extracts the isotopic patterns of the specified compound, 
#' and then performs an isotopic pattern analysis to estimate the isotopic abundance of the labelling isotope. 
#' This is performed through a weighted non-linear least squares fitting procedure, 
#' where the resulting estimate is the value for which the theoretical pattern best reproduces 
#' the experimental one. 
#' During the fitting, the experimental signals are given weights proportional to the square root 
#' of their intensity, to correct for the non unifor variance at different intensity levels. 
#' The theoretical patterns are computed using the \code{\link[ecipex]{ecipex}} R package.
#' 
#' @section Block diagram:
#' The isotopic pattern analysis can be divided into the following steps:
#' \enumerate{
#' \item Starting from a class \code{xcmsSet} object (from the \code{xcms} R package), 
#' generate a data frame of peak signal intensities or areas, 
#' with each column corresponding to a sample. 
#' This step can be avoided if the data frame is already available (obtained by other means);
#' \item Extract from the data frame the experimental isotopic patterns of the specified compound 
#' (one pattern for each sample). 
#' In the chemical formula of the compound, the element whose abundance is unknown is called "X";
#' \item Normalise the patterns and estimate the abundance of the label 
#' through a weighted non-linear least squares fitting procedure.
#' \item Summarize the results.
#' }
#'
#' @docType package
#' @name IsotopicLabelling
#' 
#' @author Ruggero Ferrazza, Pietro Franceschi
#' 
#' @examples 
#' data(xcms_obj)
#' peak_table <- table_xcms(xcms_obj)
#' fitted_abundances <- main_labelling(peak_table, compound="X40H77NO8P", labelling="C", 
#'                                     mass_shift=0.05, RT=285, RT_shift=20, 
#'                                     chrom_width=7, initial_abundance=NA)
#' summary(fitted_abundances)
#' plot(fitted_abundances, type="patterns", saveplots=FALSE)
#' plot(fitted_abundances, type="residuals", saveplots=FALSE)
#' plot(fitted_abundances, type="summary", saveplots=FALSE)
#' save_labelling(fitted_abundances)
#' grouped_estimates <- group_labelling(fitted_abundances, groups=factor(c(rep("C12",4), rep("C13",4))))
#' # Other possible lipid compounds include:
#' # [PC34:1 + H]+. compound="X42H83NO8P", RT=475, chrom_width=10
#' #[TAG47:3 + NH4]+ (a minor species). compound="X50H94NO6", RT=891, chrom_width=7
#' 

library(xcms)
library(ecipex)
library(stringr)
library(gsubfn)

