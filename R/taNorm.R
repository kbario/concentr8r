#' @title Total Area Normalisation
#' @description Total area normalisation (TA) scales spectra so that they all have a total integral of one
#' @details `taNorm()` works by summing all elements of a row (spectrum) and then dividing each element by this sum
#' @param X The spectra intended to be normalised. Can either be a single spectrum in the form of a numerical array or multiple spectra in a numerical matrix with the rows being the spectra/samples and the columns being the ppm variables
#' @param noi The noise estimation for each spectra given in an array that is row matched to the given X variable (i.e., for a single spectrum, only one noi value is required, for a matrix, the number of noi values should match the number of rows in X)
#' @return This function assigns the normalised X argument (as X_ta) and the calculated dilution factors (as dilf_ta) to the global environment.
#' @seealso \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @family {Attribute-Based}
#' @examples
#' data(X, noi, ppm)
#' taNorm(X, noi)
#' par(mfrow = c(1,2))
#' plot(ppm, X[2,], xlim = c(9.5,0.25), xlab = "Chemical Shift (ppm)", ylab = "Intensity", main = "Non-Normalised Spectra", type = 'l', col = 'red')
#' points(ppm, X[1,], type = 'l', col = 'blue')
#' legend("topleft", legend = c("Spectra 1", "Spectra 2"), col = c("red", "blue"), lty = 1)
#' plot(ppm, X_ta[1,], xlim = c(9.5,0.25), xlab = "Chemical Shift (ppm)", ylab = "Intensity", main = "Normalised Spectra", type = 'l', col = 'red')
#' points(ppm, X_ta[2,], type = 'l', col = 'blue')
#' legend("topleft", legend = c("Spectra 1", "Spectra 2"), col = c("red", "blue"), lty = 1)
#' cat(dilf_ta)
#' @export

taNorm <- function(X, noi){
    if (is.null(noi)){
      stop('Please provide a noi argument')
    }
    if (is.null(dim(X))){
      if (length(noi)!=1){
        stop("Please provide a noi value that matches the given X argument")
      }
      if (is.null(length(X))){
        stop("Please provide a valid X variable. X is neither a matrix or an array")
      }
      cat('\033[0;34mCalculating Dilf... \033[0m')
      X[X<noi]=0
      Xa <- sum(X)
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mNormalising X... \033[0m')
      Xta <- X/Xa
    } else if (!is.null(dim(X))){
      if (length(noi)!=nrow(X)){
        stop("Please provide an X that has the same number of rows as noi has values")
      }
      cat('\033[0;34mCalculating Dilfs... \033[0m')
      Xs <- t(vapply(seq_len(nrow(X)), function(i){
        x <- X[i,]
        n <- noi[i]
        x[x<n]=0
        return(x)
      }))
      Xa <- unname(apply(Xs, 1, sum))
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mNormalising X... \033[0m')
      Xta <- t(vapply(seq_len(nrow(X)), function(x){
        X[x, ]/Xa[x]
      }))
      rownames(Xta) <- rownames(X)
    } else {
      stop("X cannot be normalised")
    }
    cat('\033[1;32mDone.\n\033[0m')
    assign("X_ta", Xta, envir = .GlobalEnv)
    assign("dilf_ta", Xa, envir = .GlobalEnv)
}
