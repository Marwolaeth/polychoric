# This file is used to construct the package
# It is not a part of `polychoric` package

required_packages <- c(
  'roxygen2',
  'pkgdown'
)
install.packages(required_packages, type = 'source', INSTALL_opts = '--byte-compile')

library(devtools)
package_info()
packageVersion("devtools")
packageVersion("roxygen2")

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

document()
use_rcpp_eigen()
