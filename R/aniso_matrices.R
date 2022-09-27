#' @export
det1_gamma <- function(beta) {
  1/2 * (sqrt(beta^2 + 4) - beta)
}

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
