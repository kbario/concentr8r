#' Full Width at Half-Maximum
#' @description Calculating full width at half maximum (FWHM, aka line width). This function returns the ppm difference where peak line crosses half of the peak height. It requires one signal across all spectra within ppm ranges specified in shift.
#'
#' @param X num matrix, NMR data with rows representing spectra.
#' @param ppm num array describing chemical shift positions, its length equals to nrow(X).
#' @param sh num array(2), chemical shift range of a singlet for which fwhm is calculated
#' @param sf num, operating frequency of the spectrometer (meta$SFO1)
#'
#' @return Array of line widths in ppm. To convert from ppm to Hertz (Hz), multiply values with the spectrometer frequency (column a_SF01 in meta data frame).
#' @export
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#'
#' @examples
#' read_in(path = system.file('inst/extdata',package='concentr8r'), exp_type = list(exp=c("PROF_URINE_NOESY")), n_spec='multiple')
#' lnw <- lw(X, ppm, sh = c(-0.1, 0.1), meta$a_SF01)
#' bad <- which(lnw>1)
#' if (length(bad)!=0){cat("The spectra", bad, "have line widths over 1")}else{cat("All Spectra have line widths under 1")}

lw <- function (X, ppm, sh = c(-0.1, 0.1), sf){
    idx <- get_idx(sh, ppm)
    asign = sign(diff(ppm[seq_len(2)]))
    fwhm <- apply(X[, idx], 1, function(x, pp = ppm[idx], as = asign) {
      if (min(x) < 0) {
        x <- x + abs(min(x))
        baseline = 0
      } else {
        baseline <- min(x)
      }
      height <- max(x) - baseline
      hw <- baseline + (0.5 * height)
      f <- approxfun(pp, x)
      x_new <- seq(pp[1], pp[length(pp)], by = 1e-05 * as)
      y_new <- f(x_new)
      diff(sort(abs(x_new[range(which(y_new > hw))])))
    })
    return(fwhm * sf)
}
