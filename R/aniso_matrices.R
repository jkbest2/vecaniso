#' @export
det1_gamma <- function(beta) {
  1/2 * (sqrt(beta^2 + 4) - beta)
}

#' @title Construct anisotropy matrix
#'
#' Constructs an anisotropy matrix \code{H} as a mixture of an isotropic effect
#' and an anisotropic effect. The direction of the anisotropic component is
#' determined by \code{theta}, and the its intensity by \code{beta}. The
#' isotropic mixture component is determined by \code{gamma}. Because the
#' marginal variance of the random field will vary with the determinant of
#' \code{H}, it is also possible to calculate \code{gamma} so that the
#' determinant of \code{H} is always one.
#'
#' @param theta Primary direction of anisotropy (radians)
#' @param beta Anisotropic mixture component or intensity
#' @param gamma Isotropic mixture component; optional
#' @param det1 Use value of `gamma` so that H has determinant 1?
#' @export
aniso_h <- function(theta, beta, gamma = NULL, det1 = TRUE) {
  if (is.null(gamma) && det1) {
    gamma <- det1_gamma(beta)
  }
  if (is.null(gamma)) {
    stop("Need a value for gamma")
  }
  beta >= 0 || stop("Need beta > 0")
  gamma > 0 || stop("Need gamma > 0")
  v <- c(cos(theta), sin(theta))
  diag(c(gamma, gamma)) + beta * v %*% t(v)
}
