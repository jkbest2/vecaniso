##' Standard formula interface allows easy specification of fixed effects,
##' response, and offset.
##'
##' @title Construct the data list to pass to a vecaniso model
##' @param f A formula with the response, fixed effects and offset specified
##' @param data Spatial data frame
##' @param theta vector of anisotropy directions (radians)
##' @param beta Vector of anisotropy relative intensities
##' @param mesh INLA mesh
##' @return \code{list} suitable to pass to vecaniso model
##' @author John K Best
va_data <- function(f, data, theta, beta, mesh) {
  n_obs <- nrow(data)
  n_tri <- nrow(mesh$graph$tv)
  stopifnot(length(theta) == n_tri,
            length(beta) == n_tri)

  ## Response, fixed effects, and offset are extracted from the provided data
  ## and formula
  modframe <- model.frame(f, data)
  obs <- model.response(modframe)
  if (is.null(model.offset(modframe))) {
    offs <- rep(1, n_obs)
  } else {
    offs <- model.offset(modframe)
  }
  X <- model.matrix(f, modframe)

  ## Create the spatial projection matrix
  A_spat <- va_projection(data, mesh)

  ## Use the mesh to generate the requried FEM matrices
  spde <- va_spde(theta, beta, mesh)

  list(obs = model.response(modframe),
       offset = offs,
       X = X,
       A_spat = A_spat,
       spde = spde)
}

##' Generate FEM matrices for anisotropic SPDE. Code adapted from TMB anisotropy
##' example.
##'
##' @title Generate anisotropic FEM matrices
##' @param theta vector of anisotropy directions (radians)
##' @param beta Vector of anisotropy relative intensities
##' @param mesh INLA mesh
##' @return List of data expected for \code{spde_aniso_t} TMB object
##' @author John Best
##' @export
va_spde <- function(theta, beta, mesh) {
  inla_spde <- INLA::inla.spde2.matern(mesh, alpha = 2)

  ## whi
  dset <- 1:2

  # Triangle info
  TV <- mesh$graph$tv # Triangle to vertex indexing
  V0 <- mesh$loc[TV[, 1], 1:2] # V = vertices for each triangle
  V1 <- mesh$loc[TV[, 2], 1:2]
  V2 <- mesh$loc[TV[, 3], 1:2]
  E0 <- V2 - V1 # E = edge for each triangle
  E1 <- V0 - V2
  E2 <- V1 - V0

  # Calculate Areas
  TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
  Tri_Area <- rep(NA, nrow(E0))
  for (i in 1:length(Tri_Area)) Tri_Area[i] <- TmpFn(E0[i, ], E1[i, ]) / 2 # T = area of each triangle

  G0_inv_diag <- vapply(
    seq_len(nrow(inla_spde$param.inla$M0)),
    function(i) {
      1 / inla_spde$param.inla$M0[i, i]
    }, 0.0
  )
  G0_inv <- as(diag(G0_inv_diag), "TsparseMatrix")

  list(
    n_s      = inla_spde$n.spde,
    n_tri    = nrow(TV),
    Tri_Area = Tri_Area,
    E0       = E0,
    E1       = E1,
    E2       = E2,
    TV       = TV - 1,
    G0       = inla_spde$param.inla$M0,
    ## G0_inv   = as(diag(1 / diag(inla_spde$param.inla$M0)), "TsparseMatrix")
    G0_inv   = G0_inv,
    theta = theta,
    beta = beta
  )
}

##' Generate a projection matrix to locations in \code{data}, possible
##' zeroing out some observations. Also group by year for a spatiotemporal
##' effect.
##'
##' @title Generate a projection matrix
##' @param data Spatial data frame
##' @param mesh INLA mesh to project
##' @return A sparse projection matrix
##' @author John Best
##' @export
va_projection <- function(data, mesh, zero = FALSE) {
  ## stopifnot(st_crs(data) == mesh$crs)
  ## If effect is map'd, return an all-zero projection matrix. Should save some
  ## unnecessary multiplications?
  if (zero) {
    return(va_empty_projection(data, mesh))
  }
  loc <- sf::st_coordinates(data)
  A <- INLA::inla.spde.make.A(mesh, loc)
  Matrix::drop0(A)
}

##' @describeIn va_projection Generate empty (all-zero) projection matrix of
##'   appropriate size.
##' @export
va_empty_projection <- function(data, mesh) {
  n_obs <- nrow(data)
  Matrix::Matrix(0, nrow = n_obs, ncol = mesh$n)
}
