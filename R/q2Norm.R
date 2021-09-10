#' @title Quantile/PQN Normalisation
#' @description A method of NMR spectral normalisation where the maximum intensities
#' @details The method orders spectral intensities smallest to largest and finds the mean of each 'quantile' (i.e., the mean of the maximum intensities, the mean of the 2nd maximum intensities). These means are then reassigned back to the ppm that first housed the maximum, 2nd maximum, etc., and PQN is performed on each spectrum using its quantile normalised self. The calculated dilution factor is then used to normalise the original spectrum to reduce manipulation of peak shapes.
#' @family {Reference-Based}
#' @param X A numerical matrix containing the NMR spectra to be normalised. Rows should be the spectra and columns being the chemical shift variables. Cannot normalise only as single spectrum, X must be a matrix.
#' @param uv_used The **u**ni**v**ariate measure **used** to calculate the dilution coefficient. `mode` or `median` are accepted
#' @return This function assigns the normalised X argument (as X_q2) and the calculated dilution factors (as dilf_q2) to the global environment.
#' @seealso The methods paper that describes quantile normalisation: \url{https://doi.org/10.1007/s11306-011-0350-z}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' data(X)
#' q2Norm(X, uv_used = 'mode')
#' cat(dilf_q2)
#' @export

q2Norm <- function(X, uv_used = 'mode'){
    if (is.null(uv_used)){
      stop("uv_used was left blank. Please specify which univariate variable you want to use to calculate the dilf. 'mode' and 'median' accepted.")
    }
    if (!uv_used == 'mode' & !uv_used == 'median'){
      stop("uv_used not specified correctly. Only 'mode' and 'median' are accepted.")
    }
    if (is.null(dim(X))){
      stop('X must be a matrix. q2Norm() is unable to normalise a single spectrum.')
    } else if (!is.null(dim(X))){
      cat('\033[0;34mCalculating Quantile Means... \033[0m')
      Xs <- t(apply(X, 1, sort))
      Xr <- t(apply(X, 1, rank))
      Xm <- apply(Xs, 2, mean)
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mReassigning Intensity Values... \033[0m')
      Xq <- t(vapply(seq_len(nrow(Xr)), function(i){
        Xm[Xr[i,]]
      }))
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mCalculating Quotients... \033[0m')
      q <- t(vapply(seq_len(nrow(X)), function(j){
        X[j,]/Xq[j,]
      }))
      if (uv_used == "mode"){
        cat('\033[0;34mUsing the Mode... \033[0m')
        dilf <- vapply(seq_len(nrow(q)), function(y){
          i <- q[y,]
          d <- suppressWarnings(log10(i)[!is.nan(i) & !is.infinite(i) & !is.na(i)])
          den <- suppressWarnings(stats::density(d[!is.nan(d) & !is.infinite(d) & !is.na(d)]))
          dilf <- 10^(den$x[which.max(den$y)])
          return(dilf)
        })
      }
      if (uv_used == "median"){
        cat('\033[0;34mUsing the Median... \033[0m')
        dilf <- vapply(seq_len(nrow(q)), function(z){
          i <- q[z]
          dilf <- stats::median(i)
          return(dilf)
        })
      }
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mNormalising X... \033[0m')
      Xn <- t(vapply(seq_len(nrow(X)), function(a){
        X[a,]/dilf[a]
      }))
      rownames(Xn) <- rownames(X)
      cat('\033[1;32mDone.\n\033[0m')
      assign("X_q2", Xn, envir = .GlobalEnv)
      assign("dilf_q2", dilf, envir = .GlobalEnv)
    } else {
      stop("X cannot be normalised")
    }
}
