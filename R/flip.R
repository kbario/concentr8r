#' @title Flip the Spectra
#' @description Spectra produced from experiments with small scans can often appear upside-down. This function corrects this.
#' @details NMR urine spectra that are produced from a small number of scans are often incorrectly orientated due to the water signal. This function looks at the creatinine signal (ppm = 3.05) which is in all urine spectra to see if the value is positive or negative. For negative valued spectra, they are mulitplied by -1 and flipped to the correct orientation.
#' @param X The numerical matrix containing the NMR data you wish to check for orientation. The rows must contain information of one whole spectrum and the columns contain the specific chemical shift variables.
#' @param ppm An array of chemical shift variables. ppm should be column matched to the X matrix you are normalising.
#' @param sh The concatenated ppm values that define the lower and upper bounds of the creatinine signal. Default is c(3,3.1)
#' @return Returns the original X matrix but with all values with the correct sign.
#' @author \email{kylebario1@@gmail.com}
#' @family {Data_Manipulation}
#' @examples
#' read_in(path = system.file('inst/extdata',package='concentr8r'), exp_type = list(exp=c("PROF_URINE_NOESY")), n_spec = 'multiple')
#' Xf <- flip(X, ppm)
#' plot(ppm,X[1,],type='l',xlim=c(9.5,0.25),col='red',main="Disorientated NMR spectrum",ylab="X")
#' plot(ppm,Xf[1,],type='l',xlim=c(9.5,0.25),col ='blue',main="Orientated NMR spectrum",ylab="Xf")
#' @export

flip = function(X, ppm, sh = c(3, 3.1)){
    if (is.null(dim(X))){
      if (is.null(length(X))){
        stop("Please provide an appropriate value for X. Currently, X is neither a matrix nor an array.")
      } else if (!is.null(length(X))){
        if (length(ppm)!=length(X)){
          stop("Please provide a variable for ppm and X that are equal in length")
        } else {
          id <- get_idx(sh, ppm)
          if (sum(X[id])<0){
            X <- X*-1
          }
        }
      }
    } else if (!is.null(dim(X))){
      if (length(ppm)!=ncol(X)){
        stop("Please provide a variable for ppm with the same length as X has columns")
      } else {
        id = get_idx(sh, ppm)
        X = t(apply(X, 1, function(x, idx=id){
          if (sum(x[idx])<0) {
            x=x*-1
          }
          return(x)
        }))
      }
    } else {
      stop('X cannot be flipped')
    }
    return(X)
}

