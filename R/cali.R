#' @title Chemical Shift Calibration
#' @export
#' @description Chemical shift calibration with a reference signal.
#' @param X num matrix, NMR data matrix with rows representing spectra and columns representing chemical shift position
#' @param ppm num vector, matched to columns of X
#' @param type str or num. Str: Either 'tsp' or 'glucose' for urine or blood-derived spectra, respectively (see Details). Num: ppm range of max height signal that will be used to reference to zero
#' @details Spectral calibration to a selected chemical shift reference signal.
#'    `type='tsp'`: calibration to 0 ppm using the highest peak located in interval 0 +/- 0.20 ppm (Trimethylsilylpropanoic acid resonance)
#' @return num matrix X, calibrated NMR data matrix.
#' @references Dona, A.C., \emph{et al.} (2014) Precision high-throughput proton NMR spectroscopy of human urine, serum, and plasma for large-scale metabolic phenotyping. \emph{Analytical Chemistry}. 86.19. 9887-94.
#' @examples
#' read_in(path = system.file('extdata',package='concentr8r'),
#'         exp_type = list(exp=c("PROF_URINE_NOESY")),
#'         n_spec = 'multiple')
#' plot(ppm, X[2,], type = 'l', col = 'red',
#'     main = 'NMR Spectra (Processed vs. Non Processed)',
#'     xlab = 'Chemical Shift (ppm)', ylab = 'Intensity', xlim = c(10,-1))
#' Xc=cali(X, ppm, type='tsp')
#' points(ppm, X[2,], type = 'l', col = 'blue')
#' legend('topleft', legend = c("Uncalibrated", "Calibrated"), col = c('red', 'blue'), lty = 1)
#'
#' @author \email{torben.kimhofer@@murdoch.edu.au}


cali <- function(X, ppm, type = c('tsp')){
    if(all(is.character(type)) && (any(is.na(type[1])) || length(type)>1)){stop('Check function argument `type`!')}
    if(all(is.numeric(type)) && length(type)!=2){stop('Check function argument `type`!')}
    if(type[1] == "tsp"){type=c(-0.2, 0.2)}
    rnam <- rownames(X)
    cnam <- colnames(X)
    if (is.numeric(type)[1]) {
      idx <- get_idx(type, ppm)
      zeroPpm <- which.min(abs(ppm[idx]-mean(type)))
      Xc<-t(vapply(seq_len(nrow(X)), function(i){
        iCorr <- zeroPpm - which.max(X[i, idx])
        if (iCorr<0) { x <- c(X[i, -c(seq_len(abs(iCorr)))], rep(0, abs(iCorr))) }
        if (iCorr>0) { x <- c(rep(0, abs(iCorr)), X[i,]) }
        if (iCorr==0) { x <- X[i,] }
        x[seq_len(length(ppm))]
      }, FUN.VALUE=ppm))
    }
    rownames(Xc) <- rnam
    colnames(Xc) <- cnam
    return(Xc)
}
