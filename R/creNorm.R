#' @title Creatinine Normalisation
#' @description Creatinine Normalisation (CN) is a useful method much like region of interest normalisation that can normalise spectra based on the total area of the creatinine signal at the chemical shift 3.05ppm.
#' @details `creNorm()` works by dividing each element in a row with the sum of the values from its Creatinine signal.
#' @family {Attribute-Based}
#' @param X The spectra intended to be normalised. Can either be a single spectrum in the form of a numerical array or multiple spectra in a numerical matrix with the rows being the spectra/samples and the columns being the ppm variables
#' @param ppm A numerical array holding the chemical shift values of the X matrix. Only necessary when X is an array, not when X is a matrix
#' @param cre3 A concatenated numerical value of the lower and upper ppm values where the creatinine peak at 3.05 starts and ends.
#' @param cre4 A concatenated numerical value of the lower and upper ppm values where the creatinine peak at 4.05 starts and ends.
#' @param err The level of error given when calculating the creatinine peak ratios. interperted as a percentage (i.e., 5 = 5%)
#' @return This function assigns the normalised X argument (as X_cre) and the calculated dilution factors (as dilf_cre) to the global environment.
#' @author \email{kylebario1@@gmail.com}
#' @seealso More on the methodology of CN and issue with using it are outlined here: \url{https://doi.org/10.1021/ac051632c}
#' @examples
#' # When X contains multiple spectra, ppm is not required
#' data(X, ppm)
#' creNorm(X)
#' cat(dilf_cre)
#'
#' # When X has only one spectrum, ppm is required
#' data(X, ppm)
#' creNorm(X[1,], ppm)
#' cat(dilf_cre)
#'
#' @export

creNorm <- function(X, ppm = NULL, cre3 = c(3, 3.1), cre4 = c(4, 4.1), err = 5){
    if (length(cre3)!=2 | length(cre4)!=2){
      stop("Please provide only two values for each of the args cre3 and cre4. The first for each should be the lower bounds of the creatinine regions and the second should be the upper bounds.")
    }
    if (is.null(err)){
      err = 5
    }
    if (is.null(dim(X))){
      if (is.null(length(X))){
        stop("Please provide a valid X variable. X is neither a matrix or an array")
      }
      if (is.null(ppm)){
        stop("Please provide a X-matched ppm. None was provided and ppm cannot be determined from a single spectrum")
      } else if (length(ppm)!=length(X)){
        stop('Please make sure that the length of X and ppm match')
      }
      cat('\033[0;34mCalculating Dilfs... \033[0m')
      i3 <- shift_pickr(X, ppm, cre3, 0.005)
      i4 <- shift_pickr(X, ppm, cre4, 0.005)
      a3 <- sum(X[i3])
      a4 <- sum(X[i4])
      r <- a4/a3
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mChecking creatinine peak ratio... \033[0m')
      er <- ((2/3)/100)*err
      lo <- (2/3)-er
      up <- (2/3)+er
      if(r<=up & r>=lo){
        cat('\033[1;32mRatio is within limit.\n\033[0m')
      } else{
        cat('\033[1;31mThe provided spectra is outside of the error limit.\n\033[1;33mRefer to Df_cre for more information.\n\033[0m')
      }
      e <- as.array(r<=up & r>=lo)
      df <- data.frame(a3, r, e)
      colnames(df) <- c('dilf', 'ratio', paste('ratio within a', err, '% error margin'))
      cat('\033[0;34mNormalising X... \033[0m')
      Xn <- X/a3
    } else if (!is.null(dim(X))){
      cat('\033[0;34mCalculating Dilfs... \033[0m')
      if (is.null(ppm)){
        p <- as.numeric(colnames(X))
      } else {
        if (length(ppm)!=ncol(X)){
          stop('Please provide a column-matched ppm and X variable')
        } else {
          p <- ppm
        }
      }
      i3 <- shift_pickr(X, p, cre3, 0.005)
      i4 <- shift_pickr(X, p, cre4, 0.005)
      a3 <- vapply(seq_len(nrow(X)),function(i){
        j <- i3[i,]
        a <- sum(X[i,j])
        return(a)
      }, FUN.VALUE = 1.1)
      a4 <- vapply(seq_len(nrow(X)),function(i){
        j <- i4[i,]
        a <- sum(X[i,j])
        return(a)
      }, FUN.VALUE = 1.1)
      r <- a3/a4
      cat('\033[1;32mDone.\n\033[0m')
      cat('\033[0;34mChecking creatinine peak ratios... \033[0m')
      er <- ((2/3)/100)*err
      lo <- (2/3)-er
      up <- (2/3)+er
      if(all(r<=up & r>=lo)){
        cat('\033[1;32mAll within limits.\n\033[0m')
      } else{
        cat('\033[1;31mspec', which(r>=up | r<=lo), 'are outside error limits.\n\033[1;33mRefer to Df_cre for more information.\n\033[0m')
      }
      e <- as.array(r<=up & r>=lo)
      df <- data.frame(a3, r, e)
      colnames(df) <- c('dilf', 'ratio', paste('ratio within a', err, '% error margin'))
      cat('\033[0;34mNormalising X... \033[0m')
      Xn <- t(vapply(seq_len(nrow(X)), function(i){
        X[i,]/a3[i]
      }, FUN.VALUE = X[1,]))
      rownames(Xn) <- rownames(X)
    } else {
      stop("X cannot be normalised")
    }
    cat('\033[1;32mDone.\n\033[0m')
    assign("X_cre", Xn, envir = .GlobalEnv)
    assign("dilf_cre", a3, envir = .GlobalEnv)
    assign("Df_cre", df, envir = .GlobalEnv)
}
