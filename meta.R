# This file is used to construct the package
# It is not a part of `polychoric` package

required_packages <- c(
  'roxygen2',
  'pkgdown'
)
install.packages(required_packages, type = 'source', INSTALL_opts = '--byte-compile')
library(devtools)

use_git()
use_github()

load_all()

?check
check(document = TRUE, remote = TRUE)

use_mit_license()

document()

use_data(gss12_values)
use_data_raw('gss_extract_values')
