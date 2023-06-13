##' Use the data, parameters, map, and random objects to construct an
##' autodifferentiated objective function that can be optimized.
##'
##' @title Construct objective function
##' @param data Data list from \code{\link{va_data}}
##' @param pars Parameter list from \code{\link{va_pars}}
##' @param map Map list from \code{\link{va_map}}
##' @param random Random parameter vector from \code{\link{va_random}}
##' @return A TMB ADFun suitable for fitting
##' @author John K Best
##' @export
va_obj <- function(data, pars, map = list(), random = c("omega")) {
  TMB::MakeADFun(data = data,
                 parameters = pars,
                 map = map,
                 random = random,
                 DLL = "vecaniso")
}

##' @describeIn va_obj Simulate from the model
##' @export
va_sim <- function(obj) {
  obj$simulate()
}

##' @describeIn va_obj Fit the model using \code{nlminb}
##' @export
va_fit <- function(obj) {
  nlminb(obj$par, obj$fn, obj$gr)
}

##' @describeIn va_obj Get the REPORTed values
##' @export
va_report <- function(obj) {
  obj$report()
}

##' @describeIn va_obj Get the \code{\link[TMB]{sdreport}}
##' @export
va_sdreport <- function(obj, ...) {
  TMB::sdreport(obj, ...)
}
