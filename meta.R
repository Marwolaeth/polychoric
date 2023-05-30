required_packages <- c(
  'roxygen2',
  'pkgdown'
)
install.packages(required_packages, type = 'source', INSTALL_opts = '--byte-compile')

library(devtools)
package_info()
packageVersion("devtools")
packageVersion("roxygen2")
# library(RcppEigen)
# 
use_rcpp_eigen()
use_package("RcppEigen")
document()
# Rcpp::compileAttributes()
# 
#
use_git()
use_github()
# 
# usethis::use_testthat()
# use_test()

load_all()
# .rs.setClangDiagnostics(0)

?check
check(document = TRUE, remote = TRUE)

use_mit_license()
