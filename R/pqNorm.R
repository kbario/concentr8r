#' @title Probabilistic Quotient Normalisation
#' @description PQN is currently the gold standard method used to normalise NMR spectra.
#' @details `pqNorm()` functions by calculating the median value for each column and divides each element of a spectrum by its corresponding median. The most frequently occurring quotient (result of division) is then called the dilution factor for that spectrum and all elements within that spectrum are divided by it, scaling the spectrum.
#' @family {Reference-Based}
#' @param X The numerical matrix containing the NMR data you wish to normalise. **This should be a preprocessed matrix** with baseline correction, tsp calibration and non-quantitative region removal performed on it. The rows must contain information of one whole spectrum and the columns contain the specific chemical shift variables.
#' @param noi Takes an array that is row matched to the X matrix you are normalising with the values equaling the maximum noise estimation for each spectra respectively.
#' @param use_ta Requires a boolean `TRUE` or `FALSE` if total area normalisation should be performed on the spectra before PQN is.
#' @param uv_used PQN utilises finding the median or the mode, which are both *U*ni*v*ariate methods. Recognises either the string 'median' or 'mode' to instruct which method to use. Default = 'mode'
#' @param calc_region The lower and upper bounds of the spectrum that will be used to calculate the dilution coeficient
#' @param bin_width The width of the bin when the spectra are binned
#' @return This function assigns the normalised X argument (as X_pqn) and the calculated dilution factors (as dilf_pqn) to the global environment.
#' @seealso The methods paper first describing PQN can be found here: \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' data(X, noi)
#' pqNorm(X, noi)
#' cat(dilf_pqn)
#' @export

pqNorm <- function(X, noi, use_ta = FALSE, uv_used = 'mode', calc_region = c(0.5,9.5), bin_width = 0.01){
    if (is.null(dim(X))){
      stop('X must be a matrix. pqNorm() is unable to normalise a single spectrum.')
    } else if (!is.null(dim(X))){
      if (is.null(use_ta)){
        stop("use_ta was left blank. Please specify if you would like to normalise the spectra with Total Area prior to PQN with TRUE or FALSE.")
      } else if (!is.logical(use_ta)){
        stop("Please provide a logical (TRUE or FALSE) value for use_ta")
      }
      if (is.null(uv_used)){
        stop("uv_used was left blank. Please specify which univariate variable you want to use to calculate the dilf. 'mode' and 'median' accepted.")
      } else if (!uv_used == 'mode' & !uv_used == 'median'){
        stop("uv_used not specified correctly. Only 'mode' and 'median' are accepted.")
      }
      if (length(calc_region)!=2){
        stop("Please provide only two values for calc_region. The first should be the lower bounds of the calc region and the second should be the upper bound.")
      }
      if (length(bin_width)!=1 | !is.numeric(bin_width)){
        stop("Please provide only one numerical value for bin_width")
      }
      if (!names(noi)[1]==rownames(X)[1] | !all(!is.na(match(names(noi), rownames(X)))) | !all(diff(match(names(noi), rownames(X)))==1)){
        if (length(noi) < nrow(X)){
          stop("noi does not contain as many values as X has spectra. To continue use uNorm::noise() to calculate noise estimations for your X.")
        }
        cat("\033[1;33mThe provided X and noi arguements do not match. Attempting to match them now... \033[0m")
        ma <- match(rownames(X),names(noi))
        noi_matched <- noi[ma]
        if (all(names(noi_matched)==rownames(X))){
          assign('noi_matched', noi_matched, envir = .GlobalEnv)
          cat('\033[1;32mDone.\n\033[0;34mFind noi_matched in the global environment.\033[0m')
        }
        if (!all(names(noi_matched)==rownames(X))){
          cat('\n\033[0;31mX and noi could not be matched. Please match them before applying PQN.\n\033[0m')
        }
      }
      if (use_ta){
        cat('\033[0;34mPerforming Total Area Normalisation... \033[0m')
        X <- t(vapply(seq_len(nrow(X)), function(x){
          (X[x,])/(sum(X[x,]))
        }, FUN.VALUE = X[1,]))
        cat('\033[1;32mDone.\n\033[0m')
      }
      p <- as.numeric(colnames(X))
      cat('\033[0;34mPreparing Spectra and Reference...\n\033[0m')
      idx <- get_idx(c(calc_region[1],calc_region[2]), p)
      cat('\033[0;34mSelecting ppm and Removing Noise... \033[0m')
      Xs <- t(vapply(seq_len(nrow(X)), function(i){
        n <- noi[i]
        Xs <- X[i,idx]
        Xs[Xs<n]=0
        return(Xs)
      }, FUN.VALUE = X[1,idx]))
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mBinning... \033[0m')
      Xsb <- binin(Xs, p[idx], width = bin_width)
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mCalculating Reference Spectrum... \033[0m')
      Xm <- apply(Xsb, 2, stats::median)
      cat('\033[1;32mDone.\n')
      cat('\033[0;34mCalculating Dilfs... \033[0m')
      q <- t(vapply(seq_len(nrow(Xsb)), function(j){
        Xsb[j,]/Xm
      }, FUN.VALUE = Xsb[1,]))
      if (uv_used == "mode"){
        cat('\033[0;34mUsing the Mode... \033[0m')
        dilf <- vapply(seq_len(nrow(q)), function(y){
          i <- q[y,]
          d <- suppressWarnings(log10(i)[!is.nan(i) & !is.infinite(i) & !is.na(i)])
          den <- suppressWarnings(stats::density(d[!is.nan(d) & !is.infinite(d) & !is.na(d)]))
          dilf <- 10^(den$x[which.max(den$y)])
          return(dilf)
        }, FUN.VALUE = 1.1)
      }
      if (uv_used == "median"){
        cat('\033[0;34mUsing the Median... \033[0m')
        dilf <- vapply(seq_len(nrow(q)), function(z){
          i <- q[z,]
          dilf <- stats::median(i, na.rm = TRUE)
          return(dilf)
        }, FUN.VALUE = 1.1)
      }
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mNormalising X... \033[0m')
      Xn <- t(vapply(seq_len(nrow(X)), function(a){
        X[a,]/dilf[a]
      }, FUN.VALUE = X[1,]))
      rownames(Xn) <- rownames(X)
      cat('\033[1;32mDone.\n\033[0m')
      assign("X_pqn", Xn, envir = .GlobalEnv)
      assign("dilf_pqn", dilf, envir = .GlobalEnv)
      assign('X_pqn_median', Xm, envir = .GlobalEnv)
      assign('X_pqn_binned', Xsb, envir = .GlobalEnv)
    } else {
      stop("X cannot be normalised")
    }
}
