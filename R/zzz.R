# Set progress bar option on load
.onLoad <- function(libname, pkgname) {
    options(cli.progress_show_after = 0)
}
