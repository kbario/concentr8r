#' @title Spectral data binning
#' @export
#' @param X num matrix, NMR data with spectra in rows
#' @param ppm num array, chemical shift positions, length matches to columns in `X`
#' @param width num, bin size in ppm or NULL in case `npoints` is specified
#' @param npoints num, desired number of bins per spectrum or NULL in case  `width` is specified
#' @details Equidistant binning of spectra. Specify either `width` or `npoints` argument - if both are provided, `npoints` is used. Input argument `ppm` can be omitted if chemical shift information is encoded in the column names of the NMR matrix `X`.
#' @return Numeric matrix with spectra in rows and chemical shift variables in columns.
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @examples
#' data(X, ppm)
#' Xb <- binin(X, ppm, width = 0.01)

binin <- function (X, ppm, width = NULL, npoints=NULL){
    if (is.null(ppm) && (is.matrix(X) | is.data.frame(X)) && !is.null(colnames(X))) {
      ppm <- as.numeric(colnames(X))
    } else {
      if (ncol(X)!=length(ppm) | (any(is.na(ppm)) | any(is.infinite(ppm))))
        stop("Non-matching dimensions X matrix and ppm vector or missing values in ppm.")
    }
    if (is.vector(X)) {
      X <- t(X)
    }
    if (!is.null(width) & !is.null(npoints)) {
      stop("Please specify one or the other: bin width or desired number of data points per spectrum.\n")
    }
    if (is.null(width) & is.null(npoints)) {
      stop("Define bin width in ppm or desired number of bins.")
    }
    if (!is.null(width) & is.null(npoints)) {
      if (width <= abs(diff(ppm[seq_len(2)]))) {
        stop("Bin width equals or is smaller than the difference of neighbouring ppm points.")
      }
      res=(ppm[1] - ppm[2])
      new_res=width/round(width/res)
      step=round(width/res)
      ppm_new=seq(max(ppm), min(ppm), by=-new_res)
      iid=floor(length(ppm_new)/step)
      ybin=rep(seq(iid), each=step)
      Xb <- t(apply(X, 1, function(x, ppmt = ppm_new, ppm_fres = ppm, yb=ybin) {
        sInter <- stats::approxfun(ppm_fres, x)
        s=sInter(ppmt)
        out=sapply(seq(max(yb)), function(i){
          iidx=which(yb==i)
          sum(s[iidx])
        })
        return(out)
      }))
      ppm_bin=sapply(seq(max(ybin)), function(i){
        iidx=which(ybin==i)
        mean(ppm_new[iidx])
      })
      colnames(Xb) <- ppm_bin
      rownames(Xb) <- rownames(X)
      return(Xb)
    }
    if (!is.null(npoints) & is.null(width)) {
      if (npoints >= length(ppm)) {
        stop("Input variable npoints cannot be larger or equal than length of ppm vector.")
      }
      ppm_bin <- seq(max(ppm), min(ppm), length.out = npoints)
      iid=floor(length(ppm)/npoints)
      ybin=rep(seq(npoints), each=iid)
      Xb <- t(apply(X, 1, function(s, yb=ybin) {
        out=sapply(seq(max(yb)), function(i){
          iidx=which(yb==i)
          sum(s[iidx])
        })
        return(out)
      }))
      ppm_bin=sapply(seq(max(ybin)), function(i){
        iidx=which(ybin==i)
        mean(ppm[iidx])
      })
      colnames(Xb) <- ppm_bin
      rownames(Xb) <- rownames(X)
      return(Xb)
    }
    return(NULL)
}
