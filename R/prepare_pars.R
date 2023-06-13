##' Model parameter construction, make sure dimensions are correct
##'
##' @title Construct parameter list for vecaniso model
##' @param data Data list from \code{\link{va_prepare_data}}
##' @param mesh INLA mesh
##' @return List of initial parameter values
##' @author John K Best
##' @export
va_pars <- function(data, mesh) {
  rho_est <- min(diff(range(mesh$loc[, 1])),
                   diff(range(mesh$loc[, 2]))) / 10
  list(beta = rep(0, ncol(data$X)),
       omega = rep(0, mesh$n),
       log_kappa = log(pars_kappa(rho_est)),
       log_tau = log(pars_tau(1.0, rho_est)),
       log_beta_scale = log(1.0),
       obs_lik_pars = c(1.0, 1.5)
       )
}

##' The SPDE construction uses parameters kappa and tau, while the usual Matern
##' parameters are sigma and rho. The kappa parameter is \code{sqrt(8) / rho}.
##'
##' @describeIn pars_tau Calculate kappa from rho
##' @export
pars_kappa <- function(rho) {
  sqrt(8) / rho
}

##' The SPDE construction uses parameters kappa and tau, while the usual Matern
##' parameters are sigma and rho. These functions allow converting between the
##' parameterizations.
##'
##' The tau parameter rescales the precision matrix to a given marginal standard
##' deviation, while \code{kappa} controls the correlation decay rate. In particular,
##'
##' - \code{tau = 2 * sqrt(pi) * sig * kappa},
##' - \code{kappa = sqrt(8) / rho},
##' - \code{sig = tau / (2 * sqrt(pi) * kappa)}, and
##' - \code{rho = sqrt(8) / kappa}.
##'
##' All assume that the Matern smoothness parameter is 1 and the domain is
##' two-dimensional
##'
##' @title Calculate tau from rho and sigma^2
##' @param sig Marginal standard deviation
##' @param rho Correlation range parameter
##' @param kappa Correlation decay parameter
##' @param tau Ratio of precision and correlation range
##' @return Parameter value
##'
##' @author John Best
##' @export
pars_tau <- function(sig, rho) {
  2 * sig * sqrt(pi) * pars_kappa(rho)
}

##' @describeIn pars_tau Calculate the correlation range
##' @export
pars_rho <- function(kappa) {
  sqrt(8) / kappa
}

##' @describeIn pars_tau Calculate the marginal standard deviation
##' @export
pars_sig <- function(tau, kappa) {
  tau / (2 * sqrt(pi) * kappa)
}
