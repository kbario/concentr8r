% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/noise.R
\name{noise}
\alias{noise}
\title{Maximum Noise Estimation}
\usage{
noise(X_OG, ppm_OG, sh = c(9.5, 11), sd_mult = 5, method = "simple")
}
\arguments{
\item{X_OG}{The numerical matrix containing the NMR data you wish to estimate the noise of. The rows must contain information of one whole spectrum and the columns contain the specific chemical shift variables.}

\item{ppm_OG}{An array of chemical shift variables. ppm should be column matched to the X matrix you using.}

\item{sh}{The concatenated ppm values that define the lower and upper bounds of the noise region. Default is c(9.5,11)}

\item{sd_mult}{The value that the standard deviation will be multiplied by. Default = 5. This does not need to be changed.}

\item{method}{Takes the strings 'simple' or 'topspin'. Defines what method is used to calculate the noise. Simple is the method detailed in details. Topspin is the method used by topspin which can be found in the help section of topspin.}
}
\value{
This function returns an estimation of the maximum noise level for each spectra.
}
\description{
This function estimates the maximum level of noise for each spectra in an X matrix. The output of this function can be used as the noi variable in both \code{\link[=pqNorm]{pqNorm()}} and \code{\link[=hmNorm]{hmNorm()}}
}
\details{
This function uses a combination of the functions mean and sd to find the standard deviation of the noise region as well as it's mean. The standard deviation is multiplied by five to capture any extreme outliers and added to the mean noise level to produce an estimation of maximum noise. The output of this function is intended for use in both the \code{\link[=pqNorm]{pqNorm()}} and \code{\link[=hmNorm]{hmNorm()}} functions so they can remove the noise from the spectra and produce clear results.
}
\examples{
read_in(path = system.file('extdata',package='concentr8r'),
        exp_type = list(exp=c("PROF_URINE_NOESY")),
        n_spec = 'multiple')
noi <- noise(X, ppm, c(9.5,11))
cat(noi)
}
\seealso{
Simple noise estimation methodology was adapted from histogram matching methods paper which can be found here: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}

Other {estimation}: 
\code{\link{get_idx}()},
\code{\link{lw}()},
\code{\link{shift_pickr}()}
}
\author{
\email{kylebario1@gmail.com}
}
\concept{{estimation}}
