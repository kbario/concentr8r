#' @title Region of Interest Normalisation
#' @description This function normalises the spectra based on a specific area of the spectra. This is helpful when one region is constant across all spectra.
#' @details `roiNorm()` operates similar to `creNorm()` except the normalising factor is not specificly creatinine, it can be what ever region is of interest.
#' @family {Attribute-Based}
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @param ppm A numerical array holding the chemical shift values of the X matrix. Only necessary when X is an array, not when X is a matrix
#' @param sh The numerical values defining the lower and upper regions of the Region of Interest. default = c(3,3.1).
#' @return This function assigns the normalised X argument (as X_roi) and the calculated dilution factors (as dilf_roi) to the global environment.
#' @seealso A description of Region of Interest Normalisation can be found in this paper: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' data(X)
#' roiNorm(X, sh = c(1.5,2))
#' cat(dilf_roi)
#' @export

roiNorm <- function(X, ppm = NULL, sh = c(2.5,2.75)){
  if (length(sh)!=2){
    stop("Please provide only two values for sh. The first should be the lower bounds of the region of interest and the second should be the upper bound.")
  }
  if (is.null(dim(X))){
    if (is.null(length(X))){
      stop("Please provide a valid X variable. X is neither a matrix or an array")
    }
    if (is.null(ppm)){
      stop("Please provide a X-matched ppm. None was provided and ppm cannot be determined from a single spectrum")
    }
    p <- ppm
    i <- metabom8::get_idx(sh, p)
    cat('\033[0;34mCalculating Dilf... \033[0m')
    dilf <- sum(X[i])
    cat('\033[1;32mDone.\n\033[0m')
    cat('\033[0;34mNormalising X... \033[0m')
    Xn <- X/dilf
  } else if (!is.null(dim(X))){
    if (is.null(ppm)){
      p <- as.numeric(colnames(X))
    } else {
      if (length(ppm)!=ncol(X)){
        stop('Please provide a column-matched ppm and X variable')
      } else {
        p <- ppm
      }
    }
    i <- metabom8::get_idx(sh, p)
    cat('\033[0;34mCalculating Dilfs... \033[0m')
    dilf <- sapply(1:nrow(X), function(y){
      (sum(X[y,i]))
    })
    cat('\033[1;32mDone.\n\033[0m')
    cat('\033[0;34mNormalising X... \033[0m')
    Xn <- t(sapply(1:nrow(X), function(x){
      (X[x,])/(sum(X[x,i]))
    }))
    rownames(Xn) <- rownames(X)
  } else {
    stop("X cannot be normalised")
  }
  cat('\033[1;32mDone.\n\033[0m')
  assign("X_roi", Xn, envir = .GlobalEnv)
  assign("dilf_roi", dilf, envir = .GlobalEnv)
}
