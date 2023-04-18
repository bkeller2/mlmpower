# General Imports
#' @importFrom stats aggregate as.formula coefficients na.omit
#' @importFrom stats qnorm rnorm runif sd setNames sigma

# Set progress bar option on load
.onLoad <- function(libname, pkgname) {
    options(cli.progress_show_after = 0)
}
