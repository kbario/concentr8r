#' @title Streamlined 1D NMR Preprocessing
#' @description [preproc()] is a function aimed at streamlining the preprocessing stage of analysing NMR spectra. It harnesses the power of all preprocessing functions with in concentr8r to quickly and easily process spectra.
#' @details [preproc()] carries out a range of functions. Below is a simplistic run through. See the vignette "Preprocessing" to learn more about the process.
#' # The Pipeline
#' This function streamlines the preprocessing of NMR urine spectra by combining a range of functions. It:
#' 1. Orientates the spectra correctly,
#' 2. Calibrates the spectra by a specific peak,
#' 3. Calculates the line widths of the peaks and returns a warning with the spectra that exceed the specified threshold,
#' 4. Removes the lower, upper, water and urea regions of the spectra,
#' 5. Corrects the baseline of the spectra using asymmetric least squares
#' 6. Verifies that the resulting X, ppm and meta objects match appropriately.
#' @param X A matrix containing the non-preprocessed NMR spectral data. The rows should containing all values of a single experiment, and the columns, the values of the chemical shift variables.
#' @param ppm An array of the chemical shift variables, column matched to X.
#' @param meta The matrix of metadata pertaining to the X matrix. This is crucial for the TSP calibration and line width calculation.
#' @param flip Default is set to `TRUE`. **NOT CRUCIAL FOR STANDARDLY ACQUIRED SPECTRA** This function checks the orientation of spectra. This is a particularly important function for spectra derived from a **small number of scans**; the NMR orientates spectra based on what side of the y-axis has the most signal and water has a very large impact on this and without a large amount of metabolite signal on the correct side of the spectra (as seen in spectra derived from small scans) the spectra will be orientated the wrong way.
#' @param cali Default is set to `TRUE`. This calls on the `cali()` function of concentr8r which ensures peaks are aligned based on the TSP signal. Calibration is standardly performed by the NMR but this can be performed again for certainty.
#' @param calib This is the signal you wish to use as the calibrant. For urine spectra, *tsp* is the main calibrant.
#' @param linw This argument defines the maximum line width that the tsp should have which is spectral validation. It uses the lw() function to calculate the line width of peaks. Spectra with peaks that have small line widths have sharper and more precise results which is more desirable. A maximum cutoff of 1 ensures spectra contain robust results. Consider omitting spectra with line widths over 1.
#' @param bline Default is set to `TRUE`. This argument calls on `bl()`, a baseline correcting function to smooth the spectral baselines and remove the influence of broad peaks.
#' @param lowCut A single floating point number defining the ppm value that the lower limit of the spectra are trimmed to.
#' @param uppCut A single floating point number defining the ppm value that the upper limit of the spectra are trimmed to.
#' @param watCut The lower and upper ppm values concatenated, from which the water region will be trimmed and omitted. Water regions provide no important information and should be removed prior to data analysis. Default is set to `c(4.5,5)`
#' @param ureCut The lower and upper ppm values concatenated, from which the urea region will be trimmed and omitted. Urea regions also provide no important information and should be removed prior to data analysis. Default is set to `c(5.6,6)`
#' @param noi_sh The shift of the noise region used to calculate estimation of noise
#' @return This function returns a list with:
#' 1. The processed X matrix in the first element,
#' 2. The processed ppm array in the second element, and
#' 3. The line width results in a data frame in the third element.
#' * Following the example below will extract the results quickly and easily.
#' @export
#' @author \email{kylebario1@@gmail.com}
#' @family {preproc}
#' @examples
#' read_in(path=system.file('extdata',package='concentr8r'),
#'         exp_type=list(exp=c("PROF_URINE_NOESY")),
#'         n_spec = 'multiple')
#' plot(ppm, X[2,], type = 'l', col = 'red',
#'     main = 'NMR Spectra (Processed vs. Non Processed)',
#'     xlab = 'Chemical Shift (ppm)', ylab = 'Intensity', xlim = c(10,-1))
#' preproc(X, ppm, meta, flip = TRUE, cali = TRUE, calib = 'tsp')
#' points(ppm, X[2,], type = 'l', col = 'blue')
#' legend('topleft', legend = c("Unprocessed", "Processed"), col = c('red', 'blue'), lty = 1)

preproc <- function(X, ppm, meta, bline = TRUE, flip = TRUE, cali = TRUE, calib = 'tsp', linw = 1.0, lowCut = 0.25, watCut = c(4.5,5), ureCut = c(5.6,6), uppCut = 9.5, noi_sh = c(9.5, 11)){
    if (!is.logical(bline) | is.null(bline)){
      stop("Please provide a logical (TRUE or FALSE) value for the argument bline")
    } else if (!is.logical(flip) | is.null(flip)){
      stop("Please provide a logical (TRUE or FALSE) value for the argument flip")
    } else if (!is.logical(cali) | is.null(cali)){
      stop("Please provide a logical (TRUE or FALSE) value for the argument cali")
    } else if (calib!='tsp'){
      stop("Please provide an accepted value for calib. atm 'tsp' is the only accepted value")
    } else if (!is.numeric(linw) | is.null(linw) | length(linw)!=1){
      stop("Please provide a single numerical value for the argument linw")
    } else if (!is.numeric(lowCut) | is.null(lowCut) | length(lowCut)!=1){
      stop("Please provide a single numerical value for the argument lowCut being the lower most ppm value where the spectra should be cut")
    } else if (!is.numeric(watCut) | is.null(watCut) | length(watCut)!=2){
      stop("Please provide two numerical values for the argument watCut. The first being the lower bounds of the water region as ppm and the second, the upper bounds in ppm")
    } else if (!is.numeric(ureCut) | is.null(ureCut) | length(ureCut)!=2){
      stop("Please provide two numerical values for the argument ureCut. The first being the lower bounds of the urea region as ppm and the second, the upper bounds in ppm")
    } else if (!is.numeric(uppCut) | is.null(uppCut) | length(uppCut)!=1){
      stop("Please provide a single numerical value for the argument uppCut being the upper most ppm value where the spectra should be cut")
    } else if (!is.numeric(noi_sh) | is.null(noi_sh) | length(noi_sh)!=2){
      stop("Please provide two numerical values for the argument noi_sh. The first being the lower bounds of the noise region as ppm and the second, the upper bounds in ppm")
    }

    if (is.null(dim(X))){
      if (is.null(length(X))){
        stop("Please provide an appropriate X value. X is currently neither a matrix or an array")
      } else if (length(ppm)!=length(X)){
        stop("Please provide an X and ppm value with the same length")
      } else if (!is.null(dim(meta))){
        stop("Please provide an meta argument that matches X. As X is only one spectrum consider assigning the arg meta as meta = meta[i,] where i is the same row as the spectrum is in X")
      } else if (!is.null(length(meta)) | length(meta)!= 363){
        stop("Please provide an meta variable of length 363 long")
      } else {
        X_OG <- X
        ppm_OG <- ppm

        if (flip){
          cat('\033[0;34mFlipping the spectra... \033[0m')
          Xf <- flip(X, ppm)
          cat('\033[1;32mDone.\n\033[0m')
        } else {
          cat("\033[1;33mSpectra haven't been checked for orientation.\n\033[0m")
          Xf <- X
        }

        if (cali){
          cat('\033[0;34mCalibrating to ', calib,"... ", sep = '')
          Xc <- cali(Xf, ppm, type = calib)
          cat('\033[1;32mDone.\n\033[0m')
        } else {
          cat('\033[0;34mChecking that calibration isn\'t required... \033[0m')
          if (all(meta$p_SREF_mod==1)){
            Xc <- Xf
            cat('\033[1;32mThey\'re all good.\n\033[0m')
          } else {
            cat('\033[1;33mCalibration is required\n\033[0m')
            cat('\033[0;34mCalibrating to', calib,"... ")
            Xc <- cali(Xf, ppm, type = calib)
            cat('\033[1;32mDone.\n\033[0m')
          }
        }

        cat('\033[0;34mChecking line width of spectra... \033[0m')
        lwd <- lw(Xc, ppm, sh = c(-0.1,0.1), sf = meta$a_SFO1)
        lwB <- lwd<linw
        if (all(lwB)){
          cat('\033[1;32mAll spectra have linewidths under', linw, '\n\033[0m')
        } else {
          cat('\033[1;31m', length(which(lwB==FALSE)), 'spectra have line-widths over',linw,'.\n')
          cat('\033[1;33mCheck Df_ppro for more information.\n\033[0m')
        }
        DfX <- data.frame(lwd = lwd, lwB = lwB)

        cat('\033[0;34mRemoving non-quantative regions... \033[0m')
        idx_rm <- c(get_idx(sh = c(min(ppm), lowCut), ppm), get_idx(sh = watCut, ppm), get_idx(sh = ureCut, ppm))
        # remove these indexes
        Xr <- Xc[-idx_rm]
        ppm <- ppm[-idx_rm]
        cat('\033[1;32mDone.\n\033[0m')

        if (bline){
          cat('\033[0;34mPerforming baseline correction... \033[0m')
          X <- bl(Xr)
          cat('\033[1;32mDone.\n\033[0m')
        } else {
          cat('\033[1;33mNo baseline performed.\n\033[0m')
          X <- Xr
        }

        cat('\033[0;34mCalculating Noise Estimations... \033[0m')
        noi <- noise(X, ppm, sh = c(noi_sh[1], noi_sh[2]), sd_mult = 5)
        cat('\033[1;32mDone.\n\033[0m')

        X <- X[-get_idx(sh = c(uppCut, max(ppm)), ppm)]
        ppm <- ppm[-get_idx(sh = c(uppCut, max(ppm)), ppm)]

        cat('\033[0;34mChecking that X and meta rows match... \033[0m')
        if (dim(meta)[1]==dim(X)[1]){
          cat('\033[1;32mDone.\n\033[0m')
        } else {
          stop('The number of spectra in X and meta do not match.\n')
        }
        #check X and ppm have same columns
        cat('\033[0;34mChecking that ppm length and X columns match... \033[0m')
        if (dim(X)[2]==length(ppm)){
          cat('\033[1;32mDone.\n\033[0m')
        } else {
          stop('X columns and ppm do not match.\n')
        }
      }
    } else if (!is.null(dim(X))){
      X_OG <- X
      ppm_OG <- ppm

      if (flip){
        cat('\033[0;34mFlipping the spectra... \033[0m')
        Xf <- flip(X, ppm)
        cat('\033[1;32mDone.\n\033[0m')
      } else {
        cat("\033[1;33mSpectra haven't been checked for orientation.\n\033[0m")
        Xf <- X
      }

      if (cali){
        cat('\033[0;34mCalibrating to ', calib,"... ", sep = '')
        Xc <- cali(Xf, ppm, type = calib)
        cat('\033[1;32mDone.\n\033[0m')
      } else {
        cat('\033[0;34mChecking that calibration isn\'t required... \033[0m')
        if (all(meta$p_SREF_mod==1)){
          Xc <- Xf
          cat('\033[1;32mThey\'re all good.\n\033[0m')
        } else {
          cat('\033[1;33mCalibration is required\n\033[0m')
          cat('\033[0;34mCalibrating to', calib,"... ")
          Xc <- cali(Xf, ppm, type = calib)
          cat('\033[1;32mDone.\n\033[0m')
        }
      }

      cat('\033[0;34mChecking line width of spectra... \033[0m')
      lwd <- lw(Xc, ppm, sh = c(-0.1,0.1), sf = meta$a_SFO1)
      lwB <- lwd<linw
      if (all(lwB)){
        cat('\033[1;32mAll spectra have linewidths under', linw, '\n\033[0m')
      } else {
        cat('\033[1;31m', length(which(lwB==FALSE)), 'spectra have line-widths over',linw,'.\n')
        cat('\033[1;33mCheck Df_ppro for more information.\n\033[0m')
      }
      DfX <- data.frame(lwd = lwd, lwB = lwB)
      DfX$names <- seq_len(length(rownames(X)))
      DfX$orig_names <- rownames(X)

      cat('\033[0;34mRemoving non-quantative regions... \033[0m')
      idx_rm <- c(get_idx(sh = c(min(ppm), lowCut), ppm), get_idx(sh = watCut, ppm),get_idx(sh = ureCut, ppm))
      # remove these indexes
      Xr <- Xc[,-idx_rm]
      ppm <- ppm[-idx_rm]
      cat('\033[1;32mDone.\n\033[0m')

      if (bline){
        cat('\033[0;34mPerforming baseline correction... \033[0m')
        X <- bl(Xr)
        cat('\033[1;32mDone.\n\033[0m')
      }else{
        cat('\033[1;33mNo baseline performed.\n\033[0m')
        X <- Xr
      }

      cat('\033[0;34mCalculating Noise Estimations... \033[0m')
      noi <- noise(X, ppm, sh = c(noi_sh[1], noi_sh[2]), sd_mult = 5)
      cat('\033[1;32mDone.\n\033[0m')

      X <- X[,-get_idx(sh = c(uppCut, max(ppm)), ppm)]
      ppm <- ppm[-get_idx(sh = c(uppCut, max(ppm)), ppm)]

      cat('\033[0;34mChecking that X and meta rows match... \033[0m')
      if (dim(meta)[1]==dim(X)[1]){
        cat('\033[1;32mDone.\n\033[0m')
      } else {
        stop('The number of spectra in X and meta do not match.\n')
      }
      cat('\033[0;34mChecking that ppm length and X columns match... \033[0m')
      if (dim(X)[2]==length(ppm)){
        cat('\033[1;32mDone.\n\033[0m')
      } else {
        stop('X columns and ppm do not match.\n')
      }
    } else {
      stop("X cannot be preprocessed")
    }
    assign("X", X, envir = .GlobalEnv)
    assign("ppm", ppm, envir = .GlobalEnv)
    assign("Df_ppro", DfX, envir = .GlobalEnv)
    assign("noi", noi, envir = .GlobalEnv)
    assign("X_OG", X_OG, envir = .GlobalEnv)
    assign("ppm_OG", ppm_OG, envir = .GlobalEnv)
    }
