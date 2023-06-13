##' Specified random effects to be marginalized out using Laplace approximation
##' based on whether they are estimated or not.
##'
##' @title Create random vector
##' @param map Map list from \code{\link{va_map}}
##' @return Vector of parameter names to be marginalized out
##' @author John K Best
##' @export
va_random <- function(map) {
  rand <- c("omega")
  setdiff(rand, names(map))
}
