library(devtools)

use_r('preproc')
use_r('taNorm')
use_r('creNorm')
use_r('roiNorm')
use_r('q1Norm')
use_r('q2Norm')
use_r('pqNorm')
use_r('hmNorm')
use_r('vecNorm')
use_r('xfNorm')

use_r('noise')
use_r('shift_pickr')
use_r('flip')

check()
use_mit_license("Kyle Bario")
document()
use_package('metabom8')
use_readme_rmd()

tools::checkRdaFiles('/Users/kylebario/Documents/coding/concentr8r/data/noi.rda')
tools::resaveRdaFiles('/Users/kylebario/Documents/coding/concentr8r/data/noi.rda', compress = c('bzip2'))

tools::checkRdaFiles('/Users/kylebario/Documents/coding/concentr8r/data/ppm.rda')
tools::resaveRdaFiles('/Users/kylebario/Documents/coding/concentr8r/data/ppm.rda', compress = c('bzip2'))

tools::checkRdaFiles('/Users/kylebario/Documents/coding/concentr8r/data/X.rda')
tools::resaveRdaFiles('/Users/kylebario/Documents/coding/concentr8r/data/X.rda', compress = c('bzip2'))

data("X", "ppm", "noi")

suppressWarnings(metabom8::read1d_proc(path=system.file('extdata',package='concentr8r'),exp_type=list(exp=c("PROF_URINE_NOESY"))))
preproc(X, ppm, meta, flip = T, cali = T, calib = 'tsp')
matspec(X, ppm)
use_data(X, ppm, noi, meta, osmo, overwrite = T, compress = "bzip2")


library(pkgdown)
pkgdown::build_site()
build_readme()


