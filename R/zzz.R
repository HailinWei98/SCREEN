screen_default_options <- list(
  SCREEN.Gviz = TRUE
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(x = screen_default_options) %in% names(x = op))
  if (any(toset)) options(screen_default_options[toset])
  invisible(x = NULL)
}