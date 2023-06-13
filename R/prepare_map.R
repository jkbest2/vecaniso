##' Simple mechanism for mapping (not estimating) certain parameters. Currently
##' only allows turning off whole parameter vectors.
##'
##' @title Prepare the map argument
##' @param mappars Vector of parameter names to be mapped. can be left empty
##' @param pars Parameter object from \code{\link{va_pars}}
##' @return List of mapped parameters, possibly empty
##' @author John K Best
va_map <- function(mappars = NULL, pars) {
  maplist <- list()
  for (par in mappars) {
    npar <- length(pars[[par]])
    maplist[par] <- factor(rep(NA, npar))
  }
  maplist
}
