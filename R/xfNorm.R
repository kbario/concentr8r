#' @title Normalisation based on an External Factor
#' @description This function utilises an array of dilution coefficients acquired externally to normalise the spectra.
#' @details `xfNorm()` works by dividing all values within a spectrum by the value of its external factor.
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @param xfactor A numerical array of values that correspond to the dilution coefficients of the spectra. Must be row-matched to the provided X matrix.
#' @return This function returns a list with:
#' 1. The normalised X matrix in the first list element, and
#' 2. A numerical array of the corresponding dilution factors in the second.
#' * Following the example below will extract the results quickly and easily.
#' @seealso \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @family {Reference-Based}
#' @examples
#' data(X, osmo)
#' xfNorm(X, osmo)
#' @export

xfNorm <- function(X, xfactor){
  if (is.null(dim(X))){
    if (is.null(length(X))){
      stop("Please provide a valid X variable. X is neither a matrix or an array")
    } else if (length(xfactor)!=1){
      stop("Please provide the same number of xfactors as spectra")
    } else {
      cat('\033[0;34mNormalising X... \033[0m')
      Xn <- X/xfactor
    }
  } else if (!is.null(dim(X))){
    if (length(xfactor)!=nrow(X)){
      stop("Please provide a xfactor of same length as X has spectra (rows)")
    } else {
      cat('\033[0;34mNormalising X... \033[0m')
      Xn <- t(sapply(1:nrow(X), function(x){
        X[x, ]/xfactor[x]
      }))
      rownames(Xn) <- rownames(X)
    }
  }
  cat('\033[1;32mDone.\n\033[0m')
  assign("X_xf", Xn, envir = .GlobalEnv)
  assign("dilf_xf", xfactor, envir = .GlobalEnv)
}
