\name{xcms_obj}
\alias{xcms_obj}
\docType{data}

\title{
Example of object of the class \code{xcmsSet}, from the \code{xcms} R package
}

\description{
This data object is an example of the required experimental data, that can be processed by the \code{\link{table_xcms}} function of the \code{\link{IsotopicLabeling}} R package.
 
It contains LC-MS data relative to \eqn{^{13}}C labeling experiments on 8 samples: the first 4 represent the lipid extracts of cell coltures grown under normal conditions (natural \eqn{^{13}}C abundance), whereas in the last 4 samples the lipids were extracted from cells grown using with uniformly-labeled \eqn{^{13}}C glucose (99\% \eqn{^{13}}C labeling). The experimental patterns of most lipids are consequently affected, and the \code{\link{IsotopicLabeling}} package is capable of analyzing such patterns. 

}

\usage{data("xcms_obj")}