test_that("H anisotropy matrices constructed correctly", {
  beta <- seq(0, 10, length.out = 21)

  ## Make sure that we're always getting positive gamma values
  expect_true(all(det1_gamma(beta) > 0))
  expect_true(det1_gamma(beta[1]) == 1)
  ## Check errors
  expect_error(aniso_h(pi, 1.0, det1 = FALSE))
  expect_error(aniso_h(pi, -1.0, det1 = TRUE))
  expect_error(aniso_h(pi, 1.0, -1.0, det1 = FALSE))
  ## Check that H gives determinant 1 matrices for range of beta values
  expect_true(all.equal(sapply(beta,
                               \(b) det(aniso_h(pi, b, det1 = TRUE))),
                        rep_len(1, length(beta))))
})
