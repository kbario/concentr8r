#' @title Import 1D NMR spectra
#' @export
#' @details This function imports TopSpin processed NMR spectra as well as spectrometer and processing parameters found in files \emph{acqus} and \emph{procs}. Experiments can be filtered according to data acquisition variables using the \code{exp_type} argument: For example, to read standard 1D NMR experiments use \code{exp_type=list(exp='noesygppr1d')}. More than one argument can be provided as list element.
#' @param path char, path to file directory containing spectra
#' @param exp_type named list, filter for acquisition parameters of experiments to read-in (see Details)
#' @param n_max int, maximum number of experiments to read-in
#' @param filter logic, remove experiments with incomplete file systems (TRUE is recommended)
#' @param recursive logic, if TRUE recursively search all subfolders of path for specified NMR files
#' @param verbose num, different verbose levels: 0 (no info), 1 (overview), 2 (detailed), 3 (step-by-step for debugging)
#' @param n_spec string, the number of spectra being read in. Takes either '1' or 'multiple'.
#'
#' @return
#' The function exports the following three objects into the currently active R environment (no variable assignments needed):
#' \itemize{
#'   \item X, num matrix:  NMR data, spectra in rows
#'   \item ppm, num array - chemical shift positions, length matches to columns in X
#'   \item meta, data.frame - spectrometer metadata as extracted from individual \emph{acqus} files, row-matched to X
#' }
#' Objects in the R environment with the same variable names will be overwritten.
#' @examples
#' read_in(path = system.file('inst/extdata',package='concentr8r'), exp_type = list(exp=c("PROF_URINE_NOESY")), n_spec = 'multiple')
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom stats approxfun
#' @section

read_in <- function(path, exp_type = list(exp = c("PROF_PLASMA_CPMG")), n_max = 1000, filter = TRUE, recursive = TRUE, verbose = TRUE, n_spec = "multiple") {
    if (n_spec!='multiple' & n_spec!= '1'){
      stop("Please provide a valid argument for n_spec. Must be a string, either '1' or 'multiple'.")
    }
    if (n_spec=="multiple"){
      path <- path.expand(path)
      f_list <- .detect1d_procs(path, n_max = 1e6, filter, recursive, verbose)
      pars <- .extract_pars1d(f_list)
      exp_filt <- .filterExp_files(pars, exp_type, f_list, n_max)
      f_list <- exp_filt[[1]]
      pars <- exp_filt[[2]]
      ppm_ref <- .chemShift(swidth = pars$a_SW[1], offset = pars$p_OFFSET[1], si = pars$p_SI[1])
      out <- vapply(seq_along(f_list[[1]]), function(s, ppref = ppm_ref) {
        csF2_ppm <- .chemShift(swidth = pars$a_SW[s], offset = pars$p_OFFSET[s], si = pars$p_SI[s])
        byteord <- c(little = 0, big = 1)
        spec <- readBin(f_list$f_1r[s], what = "int", n = pars$p_FTSIZE[s], size=4,
                        signed = TRUE, endian = names(byteord)[match(pars$p_BYTORDP[s], byteord)])
        spec <- (spec * (2^(pars$p_NC_proc[s])))
        nspec <- length(spec)
        f_spec <- approxfun(x = csF2_ppm, y = spec)
        spec_inter <- f_spec(ppref)
        return(spec_inter)
      }, FUN.VALUE = ppm_ref)
      out <- t(out)
      colnames(out) <- ppm_ref
      fnam <- strsplit(f_list$f_1r, .Platform$file.sep)
      idx_keep <- which((apply(do.call(rbind, fnam), 2, function(x) length(unique(x)))) > 1)
      fnam <- vapply(fnam, function(x, st = idx_keep) {
        paste(x[idx_keep], collapse = .Platform$file.sep)
      }, FUN.VALUE = "")
      rownames(out) <- fnam
      rownames(pars) <- fnam
      X[is.na(X)]=0
      assign("X", out, envir = .GlobalEnv)
      assign("ppm", ppm_ref, envir = .GlobalEnv)
      assign("meta", pars, envir = .GlobalEnv)
    }
    if (n_spec=="1"){
      readin <- function(path) {
        path = gsub(paste0(.Platform$file.sep, "$"), '', path)
        ff<-list.files(path = paste(path, 'pdata', '1', sep=.Platform$file.sep), pattern = "^1r$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
        if( length(ff) == 1 ) {
          f_list=list(list(f_procs=paste(path, 'pdata', '1', 'procs', sep=.Platform$file.sep), f_acqus=paste(path, 'acqus', sep=.Platform$file.sep)))
          meta<-extract_pars1d_(f_list)
          endianness<-ifelse(meta$p_BYTORDP!=0, 'big', 'little')
          spec <- readBin(ff, what = "int", n = meta$p_FTSIZE, size = 4, signed = T, endian = endianness)
          spec <- ((2^meta$p_NC_proc) * spec)
          swp <-  meta$p_SW_p/meta$p_SF
          dppm <- swp/(length(spec) - 1)
          offset <- meta$p_OFFSET
          ppm <- seq(from = offset, to = (offset - swp), by = -dppm)
          assign("x", spec, envir = .GlobalEnv)
          assign("p", ppm, envir = .GlobalEnv)
          assign("m", meta, envir = .GlobalEnv)
        } else {
          return(NULL)
        }
      }

    }
}

#' @title Read Bruker NMR paramter files - helper function read1d
#' @param f_list list, intact files system for NMR experiments. See fct checkFiles1d
#' @return data frame of spectrometer acquisition metadata
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @section

extract_pars1d_ <- function(f_list) {
    out <- lapply(seq(length(f_list)), function(i) {

      f_procs <- f_list[[i]]$f_procs[i]
      fhand <- file(f_procs, open = "r")
      f_procs <- readLines(fhand, n = -1, warn = FALSE)
      close(fhand)

      idx <- grep("..", f_procs, fixed = TRUE)
      f_procs[idx] <- vapply(idx, function(i) {
        gsub(" .*", f_procs[i + 1], f_procs[i])
      }, FUN.VALUE = "")
      out <- strsplit(gsub("^##\\$", "", grep("^##\\$", f_procs, value = TRUE, fixed = FALSE), fixed = FALSE), "=")
      d_procs_val <- gsub("^ ", "", vapply(out, "[[", 2, FUN.VALUE = ""))
      names(d_procs_val) <- paste0("p_", vapply(out, "[[", 1, FUN.VALUE = ""))

      f_acqu <- f_list[[i]]$f_acqus[i]
      fhand <- file(f_acqu, open = "r")
      f_acqu <- readLines(fhand, n = -1, warn = FALSE)
      close(fhand)

      idx <- grep("..", f_acqu, fixed = TRUE)
      f_acqu[idx] <- vapply(idx, function(i) {
        gsub(" .*", f_acqu[i + 1], f_acqu[i])
      }, FUN.VALUE = "")

      out <- strsplit(gsub("^##\\$", "", grep("^##\\$", f_acqu, value = TRUE, fixed = FALSE), fixed = FALSE), "=")
      d_acqu_val <- gsub("^ ", "", vapply(out, "[[", 2, FUN.VALUE = ""))
      names(d_acqu_val) <- paste0("a_", vapply(out, "[[", 1, FUN.VALUE = ""))

      idx <- grep("date", names(d_acqu_val), ignore.case = TRUE)
      d_acqu_val[idx] <- as.character(as.POSIXct(x = "01/01/1970 00:00:00", format = "%d/%m/%Y %H:%M:%S") + (as.numeric(d_acqu_val[idx])))
      pars <- c(d_acqu_val, d_procs_val)
      return(pars)
    })

    out_le <- vapply(out, length, FUN.VALUE = 1)
    if (length(unique(out_le)) > 1) {
      cnam <- unique(as.vector(vapply(out, names, FUN.VALUE = out[[1]])))
      out_df <- matrix(NA, nrow = 1, ncol = length(cnam))
      out <- as.data.frame(t(vapply(out, function(x, odf = out_df, cc = cnam) {
        odf[1, match(names(x), cnam)] <- x
        return(odf)
      }, FUN.VALUE = out[[1]])))
      colnames(out) <- cnam
    }
    if (is.list(out)) {
      out <- do.call(rbind, out)
    }
    if (nrow(out) != length(f_list)) {
      out <- t(out)
    }
    dtype_num <- apply(out, 2, function(x) {
      !any(is.na(suppressWarnings(as.numeric(x))))
    })
    out <- as.data.frame(out)
    out[, dtype_num] <- apply(out[, dtype_num], 2, as.numeric)
    return(out)
}

#' @title filtering Bruker NMR paramter files - helper function read1d
#' @param n_max int, maximum number of experiments to read-in
#' @param pars the parameters
#' @param exp_type The experiment type filtered for
#' @param f_list list, intact files system for NMR experiments. See fct checkFiles1d
#' @return data frame filtered experiment files
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @section

.filterExp_files <- function (pars, exp_type, f_list, n_max) {
    idx <- match(toupper(names(exp_type)), toupper(gsub("[ap]_", "", colnames(pars))))
    if (length(idx) == 0) {
      stop("No parameter(s) found that that match the specification. Check for exp_type argument for typos and parameter choices in acqus and procs files.")
    }
    idx_na <- which(is.na(idx))
    if (length(idx_na) > 0) {
      if (length(idx_na) == length(idx)) {
        stop("No matching paramter names found. Check input argument exp_type.")
      }
      else {
        message(paste("Experiment filter", names(exp_type)[idx_na], "not in NMR acquisition list. Using remaining arguments to filter :", names(exp_type)[-idx_na]))
        idx = idx[-idx_na]
      }
    }
    fmat <- sapply(seq(length(idx)), function(i) {
      vars = gsub("^<|>$", "", pars[, idx[i]])
      vars %in% exp_type[[i]]
    })
    idx_filt = apply(fmat == 1, 1, all)
    if (!any(idx_filt)) {
      stop("No files found that match the specified parameter specification levels.")
    }
    f_list = lapply(f_list, function(x) {
      x[idx_filt]
    })
    idx_filt = which(idx_filt)
    if (length(idx_filt) > n_max) {
      idx_filt = idx_filt[seq_len(n_max)]
      f_list = lapply(f_list, function(x, idx = seq_len(n_max)) {
        x[idx]
      })
    }
    else {
      pars <- pars[idx_filt, ]
    }
    fnam <- strsplit(f_list$f_1r, .Platform$file.sep)
    idx_keep <- which((apply(do.call(rbind, fnam), 2, function(x) length(unique(x)))) > 1)
    if (length(idx_keep) > 1) {
      rack_ <- t(vapply(seq_len(length(fnam)), function(i, iid = idx_keep[length(idx_keep) - 1]) {
        c(fnam[[i]][iid], pars$a_DATE[i])
      }, FUN.VALUE = c("", "")))
      colnames(rack_) <- c("a", "b")
      rack_order_ <- ddply(as.data.frame(rack_), .(a), function(x) {
        mean(as.POSIXct(x$b))
      })
      rord_fac <- order(rack_order_$V1) * 1e+05
      fnam1 <- vapply(fnam, function(x, st = idx_keep[length(idx_keep) -  1]) {
        x[st]
      }, FUN.VALUE = "")
      rord_fac = rord_fac[match(fnam1, rack_order_$a)]
      exp_ <- vapply(fnam, function(x, st = idx_keep[length(idx_keep)]) {
        x[st]
      }, FUN.VALUE = "")
      rr <- order(as.numeric(exp_) + rord_fac)
    }
    else {
      exp_ <- vapply(fnam, function(x, st = idx_keep) {
        x[st]
      }, FUN.VALUE = "")
      rr <- order(as.numeric(exp_))
    }
    pars = pars[rr, ]
    f_list = lapply(f_list, function(x) {
      x[rr]
    })
    pars$a_DATE = as.POSIXct(pars$a_DATE)
    return(list(f_list, pars))
}

#' @title Read Bruker NMR paramter files - helper function read1d
#' @param f_list list, intact files system for NMR experiments. See fct checkFiles1d
#' @return data frame of spectrometer acquisition metadata
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @keywords internal
#' @section

.extract_pars1d <- function(f_list) {
    out <- lapply(seq(f_list[[1]]), function(i) {
      f_procs <- f_list$f_procs[i]
      fhand <- file(f_procs, open = "r")
      f_procs <- readLines(fhand, n = -1, warn = FALSE)
      close(fhand)
      idx <- grep("..", f_procs, fixed = TRUE)
      f_procs[idx] <- vapply(idx, function(i) {
        gsub(" .*", f_procs[i + 1], f_procs[i])
      }, FUN.VALUE = "")
      out <- strsplit(gsub("^##\\$", "", grep("^##\\$", f_procs, value = TRUE, fixed = FALSE), fixed = FALSE), "=")
      d_procs_val <- gsub("^ ", "", vapply(out, "[[", 2, FUN.VALUE = ""))
      names(d_procs_val) <- paste0("p_", vapply(out, "[[", 1, FUN.VALUE = ""))

      f_acqu <- f_list$f_acqus[i]
      fhand <- file(f_acqu, open = "r")
      f_acqu <- readLines(fhand, n = -1, warn = FALSE)
      close(fhand)

      idx <- grep("..", f_acqu, fixed = TRUE)
      f_acqu[idx] <- vapply(idx, function(i) {
        gsub(" .*", f_acqu[i + 1], f_acqu[i])
      }, FUN.VALUE = "")

      out <- strsplit(gsub("^##\\$", "", grep("^##\\$", f_acqu, value = TRUE, fixed = FALSE),fixed = FALSE), "=")
      d_acqu_val <- gsub("^ ", "", vapply(out, "[[", 2, FUN.VALUE = ""))
      names(d_acqu_val) <- paste0("a_", vapply(out, "[[", 1, FUN.VALUE = ""))
      idx <- grep("date", names(d_acqu_val), ignore.case = TRUE)
      d_acqu_val[idx] <- as.character(as.POSIXct(x = "01/01/1970 00:00:00", format = "%d/%m/%Y %H:%M:%S") + (as.numeric(d_acqu_val[idx])))
      pars <- c(d_acqu_val, d_procs_val)
      return(pars)
    })
    out_le <- vapply(out, length, FUN.VALUE = 1)
    if (length(unique(out_le)) > 1) {
      cnam <- unique(as.vector(unlist(lapply(out, names))))
      out_df <- matrix(NA, nrow = 1, ncol = length(cnam))
      out <- as.data.frame(t(sapply(out, function(x, odf = out_df, cc = cnam) {
        odf[1, match(names(x), cnam)] <- x
        return(odf)
      })))
      colnames(out) <- cnam
    }
    if (!is.data.frame(out)) {
      out <- do.call(rbind, out)
    }
    if (nrow(out) != length(f_list[[1]])) {
      out <- t(out)
    }
    dtype_num <- apply(out, 2, function(x) {
      !any(is.na(suppressWarnings(as.numeric(x))))
    })
    out <- as.data.frame(out)
    out[, dtype_num] <- apply(out[, dtype_num], 2, as.numeric)
    rownames(out) <- f_list[[2]]
    return(out)
}

#' @title Check for intact file systems - helper function read1d!
#' @param datapath char, File directory containing spectra
#' @param procs_exp num, Topspin processing experiment ID
#' @param n_max int, Maximum number of spectra to read-in
#' @param filter lobic, filter for intact file systems (TRUE is recommended)
#' @return list of file descriptors: experiment folder path, folder id, procs file, acqus file, 1r file.
#' list(path=p_intact, exp_no=exp_no, f_procs=f_procs, f_acqus=f_acqus, f_1r=f_1r)
#' @author \email{torben.kimhofer@@murdoch.edu.au}
#' @keywords internal
#' @section

.detect1d_procs <- function(datapath, n_max = 10, filter = TRUE, recursive, verbose) {
    datapath <- gsub(paste0(.Platform$file.sep, "$"), "", datapath)
    f_procs <- list.files(path = datapath, pattern = "^procs$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    f_acqus <- list.files(path = datapath, pattern = "^acqus$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    f_1r <- list.files(path = datapath, pattern = "^1r$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    id_a <- gsub(paste("^", datapath, "/|/acqus$", sep = ""), "", f_acqus)
    id_p <- gsub(paste("^", datapath, "/|/pdata/.*", sep = ""), "", f_procs)
    id_f1r <- gsub(paste("^", datapath, paste0("/|/pdata/.*"), sep = ""), "", f_1r)
    idx_a <- which(id_a %in% id_f1r)
    idx_p <- which(id_p %in% id_f1r)
    idx_f1r <- which(id_f1r %in% id_a)
    if (length(idx_a) != length(id_a) || length(idx_f1r) != length(id_f1r) || length(idx_p) !=
        length(id_p)) {
      if (filter) {
        if (verbose > 1) {
          message("Reading experiments with matching acqus, procs and 1r files.")
        }
        f_acqus <- f_acqus[idx_a]
        id_a <- id_a[idx_a]
        f_procs <- f_procs[idx_p]
        id_p <- id_p[idx_p]
        f_1r <- f_1r[idx_f1r]
        id_f1r <- id_f1r[idx_f1r]
      } else {
        message("File system seems to be corrupt for some experiments. Consider function argument `filter=TRUE`")
        return(NULL)
      }
    }
    idm_f1r <- match(id_a, id_f1r)
    idm_p <- match(id_a, id_p)
    if (any(is.na(idm_f1r)) || any(is.na(idm_p))) {
      stop("check matching of this functions")
    }
    f_procs <- f_procs[idm_p]
    f_1r <- f_1r[idm_f1r]
    id_f1r <- id_f1r[idm_f1r]
    if (length(unique(c(length(f_acqus), length(f_1r), length(f_procs)))) != 1) {
      stop("Something's wrong after filtering!")
    }
    if (n_max < length(f_1r)) {
      idx_nmax <- seq_len(n_max)
      f_acqus <- f_acqus[idx_nmax]
      f_procs <- f_procs[idx_nmax]
      f_1r <- f_1r[idx_nmax]
      id_f1r <- id_f1r[idx_nmax]
      message("Reached n_max - not all spectra read-in.")
    }
    p_intact <- gsub("/acqus$", "", f_acqus)
    exp_no <- id_f1r
    return(list(path = p_intact, exp_no = exp_no, f_procs = f_procs, f_acqus = f_acqus, f_1r = f_1r))
}

#' @title Calculate chemical shift axis
#' @return num array, chemical shift for resp experiment
#' @keywords internal
#' @section
.chemShift <- function(swidth, offset, si) {
    dppm <- swidth/(si - 1)
    cshift <- seq(offset, (offset - swidth), by = -dppm)
    return(cshift)
}
