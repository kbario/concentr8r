#' @title Maximum Noise Estimation
#' @description This function estimates the maximum level of noise for each spectra in an X matrix. The output of this function can be used as the noi variable in both [pqNorm()] and [hmNorm()]
#' @details This function uses a combination of the functions mean and sd to find the standard deviation of the noise region as well as it's mean. The standard deviation is multiplied by five to capture any extreme outliers and added to the mean noise level to produce an estimation of maximum noise. The output of this function is intended for use in both the [pqNorm()] and [hmNorm()] functions so they can remove the noise from the spectra and produce clear results.
#' @param X_OG The numerical matrix containing the NMR data you wish to estimate the noise of. The rows must contain information of one whole spectrum and the columns contain the specific chemical shift variables.
#' @param ppm_OG An array of chemical shift variables. ppm should be column matched to the X matrix you using.
#' @param sh The concatenated ppm values that define the lower and upper bounds of the noise region. Default is c(9.5,11)
#' @param sd_mult The value that the standard deviation will be multiplied by. Default = 5. This does not need to be changed.
#' @param method Takes the strings 'simple' or 'topspin'. Defines what method is used to calculate the noise. Simple is the method detailed in details. Topspin is the method used by topspin which can be found in the help section of topspin.
#' @return This function returns an estimation of the maximum noise level for each spectra.
#' @seealso Simple noise estimation methodology was adapted from histogram matching methods paper which can be found here: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' read_in(path = system.file('inst/extdata',package='concentr8r'), exp_type = list(exp=c("PROF_URINE_NOESY")), n_spec = 'multiple')
#' noi <- noise(X, ppm, c(9.5,11))
#' cat(noi)
#' @family {estimation}
#' @export

noise <- function(X_OG, ppm_OG, sh = c(9.5,11), sd_mult = 5, method = 'simple'){
    if (length(sh)!=2){
      stop("Please provide only two values for sh. The first should be the lower bounds of the noise region and the second should be the upper bound.")
    }
    if (is.null(ppm_OG)){
      stop("Please provide a value for the argument ppm_OG")
    }
    if (is.null(dim(X_OG))){
      if (is.null(length(X_OG))){
        stop("Please provide a appropriate X_OG variable. The provided X_OG variable is neither a matrix nor an array")
      } else if (length(X_OG)!=length(ppm_OG)){
        stop("Please provide an X_OG and ppm_OG varibale that have the same length")
      }
    } else if (!is.null(dim(X_OG))){
      if (length(ppm_OG)!=ncol(X_OG)){
        stop("Please provide a ppm_OG variable that has same length as X_OG has columns")
      }
    }
    if (method!='simple' & method!= 'topspin'){
      stop("please provide a valid noise estimation method. 'simple' and 'topspin' are accepted")
    } else if (method=='simple'){
      if (is.null(nrow(X_OG))){
        rm <-  sd_mult*stats::sd(X_OG[get_idx(sh, ppm_OG)])+(mean((X_OG[get_idx(sh, ppm_OG)]), trim = 0.05))
      } else {
        rm <- apply(X_OG, 1, function(i){
          rm <- sd_mult*stats::sd(i[get_idx(sh, ppm_OG)])+(mean((i[get_idx(sh, ppm_OG)]), trim = 0.05))
          return(rm)
        })
      }
    } else if (method == "topspin"){
      if (is.null(nrow(X_OG))){
        x <- X_OG[get_idx(sh, ppm_OG)]
        n <- floor(length(x)/2)
        N <- n*2
        le <- length(x)
        out <- vapply(seq_len(n), function(i){
          i*(x[i]-x[le-i])
        })
        t2 <- (3*(sum(out)^2)) / (N^2 -1)
        t1 <- sum(x)^2
        rm <- sqrt((sum(x^2) - ((t1+t2)/N)) / N-1)
        return(rm)
      } else {
        rm <- t(vapply(seq_len(nrow(X_OG)), function(i){
          x <- X_OG[i, get_idx(sh, ppm_OG)]
          n <- floor(length(x)/2)
          N <- n*2
          le <- length(x)
          out <- vapply(seq_len(n), function(i){
            i*(x[i]-x[le-i])
          })
          t2 <- (3*(sum(out)^2)) / (N^2 -1)
          t1 <- sum(x)^2
          rm <- sqrt((sum(x^2) - ((t1+t2)/N)) / N-1)
          return(rm)
        }))
      }
    }
    return(rm)
}
