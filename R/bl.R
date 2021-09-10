#' @title Baseline Correction
#' @description `bl()` removes the arched baseline from an NMR spectrum making the analysis more robust
#'
#' @details NMR spectra often have broad signals produced by proteins within a sample that span a large portion of the ppm axis. These broad signals interfere with the analysis of specific peaks. By removing it, the anaylsis is more robust.
#' This is achieved by calculating the trend of the spectrum using `ptw::asysm()` and then subtracts that from the spectrum's values.
#'
#' @param x the spectrum to be baseline corrected
#'
#' @return an x with a corrected baseline
#' @importFrom ptw asysm
#' @export
#' @family {preproc}
#' @author Kyle Bario \email{kylebario1@@gmail.com}
#'
#' @examples
#' read_in(path = system.file('extdata', package = 'concentr8r'))
#' Xb <- bl_(X)

bl <- function(X){
    if (is.null(dim(X))){
      if(is.null(length(X))){
        stop("Please provide another X, this one has no dimensions or length")
      } else {
        if (any(is.na(X))){
          X[is.na(X)]=0
        }
        Xb <- X-ptw::asysm(X, maxit = 30, lambda = 1e+07)
        return(Xb)
      }
    } else if (!is.null(dim(X))){
      if (any(is.na(X))){
        X[is.na(X)]=0
      }
      Xb <- t(apply(X, 1, function(x){
        x-asysm(x, maxit = 30, lambda = 1e+07)
      }))
      return(Xb)
    }
}
