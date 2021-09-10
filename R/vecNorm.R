#' Vector Length Normalisation
#' @description This function aims to normalise NMR spectra by their vector length.
#' @details This function operates by squaring all values of the spectrum, summing them, and square rooting the result. Then each element from that spectrum is divided by this vector length.
#'
#' @param X The spectr(a/um) to be normalised. Can either be a single spectrum (in the form of an array) or multiple spectra (in the form of a matrix with columns being ppm variables and rows being the spectra)
#'
#' @return This function assigns the normalised X argument (as X_vec) and the calculated dilution factors (as dilf_vec) to the global environment.
#' @export
#'
#' @author \email{kylebario1@@gmail.com}
#' @family {Attribute-Based}
#'
#' @seealso The methods paper that describes this method of normalisation can be found here: \url{https://doi.org/10.1021/ac051632c}
#'
#' @examples
#' data(X)
#' vecNorm(X)
#' cat(dilf_vec)
#'
vecNorm <- function(X){
    if (is.null(dim(X))){
      if (is.null(length(X))){
        stop("Please provide a valid X variable. X is neither a matrix or an array")
      }
      cat('\033[0;34mCalculating Dilf... \033[0m')
      l <- sqrt(sum(X^2))
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mNormalising X... \033[0m')
      Xn <- X/l
    } else if (!is.null(dim(X))){
      cat('\033[0;34mCalculating Dilfs... \033[0m')
      l <- apply(X, 1, function(x){
        sqrt(sum(x^2))
      })
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mNormalising X... \033[0m')
      Xn <- vapply(seq_len(nrow(l)), function(y){
        X[y,]/l[y]
      })
    } else {
      stop("X cannot be normalised")
    }
    cat('\033[1;32mDone.\n\033[0m')
    assign("X_vec", Xn, envir = .GlobalEnv)
    assign("dilf_vec", l, envir = .GlobalEnv)
}
