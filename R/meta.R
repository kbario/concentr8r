#' Metadata for the NMR spectra.
#'
#' A data frame containing acquisition and processing data pertaining to the NMR spectrum.
#'
#' @format a data frame with 363 variables for two observations
#' \describe{
#'     \item{DS}{Dummy scans performed prior to acquisition}
#'     \item{NS}{Number of scans performed during the acquisition}
#'     \item{EXP}{The type of experiment/pulse sequence used to acqutire the spectrum}
#'     \item{DATE}{The date at which the sample was finished being acquired}
#'     \item{USER}{Information pertaining to the user as well as other metadata such as the success of the experiment}
#' }
"meta"
