#' Get Index
#' @description Find the elements of ppm that reside between to points
#' @details `get_idx()` asks what ppm values are between the points specified in the arg sh and returns this as an array
#'
#' @param sh The lower and upper bounds (2 numeric values) of the chemical shift region you wish to get the index for in ppm
#' @param ppm The ppm variable matched to your X
#'
#' @return A numerical array of values that can be used to index both ppm and X to specify a region to visualise or work on
#' @export
#'
#' @examples
#' # You want to visualise the Creatinine peak in your spectra
#' data(X, ppm)
#' idx <- get_idx(sh = c(3, 3.1), ppm = ppm)
#' plot(ppm[idx], X[1,idx], main = "The Creatinine Peak Found Using get_idx()", xlab = "Chemical Shift (ppm)", ylab = "Intensity", type = 'l')

get_idx <- function(sh, ppm){
    if (length(sh)!=2 | !is.numeric(sh)){
      stop("Please provide only two (2) integers in the argument sh")
    } else if (!is.null(dim(ppm)) | !is.numeric(ppm)){
      stop("Please only provide a numerical array for ppm")
    } else {
      sh <- sort(sh, decreasing = TRUE)
      which(ppm <= sh[1] & ppm >= sh[2])
  }
}
