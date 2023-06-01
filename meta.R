# This file is used to construct the package
# It is not a part of `polychoric` package
required_packages <- c(
  'roxygen2',
  'pkgdown',
  'available'
)
install.packages(required_packages, type = 'source', INSTALL_opts = '--byte-compile')
library(devtools)
available::available('polychoric') # Hey!

# use_git()
# use_github()

load_all()

check(document = TRUE, remote = TRUE, cran = TRUE)

# use_mit_license()

# use_data(gss12_values)

# use_readme_rmd()
# build_readme()

install()
