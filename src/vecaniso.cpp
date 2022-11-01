#define TMB_LIB_INIT R_init_vecaniso
#include <TMB.hpp>
#include "../inst/include/vecaniso.hpp"
#include "../inst/include/tweedielink.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // ===========================================================================
  // DATA section
  // ---------------------------------------------------------------------------
  // Observations
  DATA_VECTOR(obs);
  DATA_VECTOR(offset);

  // Fixed effects design matrices
  DATA_MATRIX(X);

  // Spatial/spatiotemporal projection matrices
  DATA_SPARSE_MATRIX(A_spat);

  // FEM matrices for SPDE spatial/spatiotemporal random effects
  DATA_STRUCT(spde, spde_vecaniso_t);

  // Observation likelihood
  // DATA_INTEGER(obs_lik);

  // ===========================================================================
  // PARAMETER section
  // ---------------------------------------------------------------------------
  // Abundance fixed effects
  PARAMETER_VECTOR(beta);
  // Abundance spatial effects
  PARAMETER_VECTOR(omega);

  // Spatial and spatiotemporal field parameters
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  // Local anisotropy scale
  PARAMETER(log_beta_scale);

  // Log catch variation parameter
  PARAMETER_VECTOR(obs_lik_pars);  // 1 for ZI log-normal, 2 for Tweedie

  // ===========================================================================
  // Derived values
  // ---------------------------------------------------------------------------
  // Number of observations
  int N_obs = obs.size();
  // Number of vertices
  int N_vert = omega.size();
  // Parameters
  Type kappa = exp(log_kappa);
  Type tau = exp(log_tau);
  Type beta_scale = exp(log_beta_scale);

  // Join negative log-likelihood accumulator. One element for each random
  // effect/random field and one for observation likelihood.
  vector<Type> jnll(2);
  jnll.setZero();

  // ===========================================================================
  // Abundance fixed effects
  // ---------------------------------------------------------------------------
  vector<Type> fixef(N_obs);
  fixef = X * beta;

  // ===========================================================================
  // Spatial random effects
  // ---------------------------------------------------------------------------
  vector<Type> spat(N_obs);
  spat = A_spat * omega;

  SparseMatrix<Type> Q_om = Q_spde(spde, kappa, beta_scale);
  SCALE_t<GMRF_t<Type>> gmrf_om = SCALE(GMRF(Q_om), tau);
  jnll(0) += gmrf_om(omega);

  SIMULATE {
    gmrf_om.simulate(omega);
    spat = A_spat * omega;

    REPORT(omega);
  }

  // ===========================================================================
  // Calculate linear predictor
  // ---------------------------------------------------------------------------
  vector<Type> linpred(N_obs);
  linpred = fixef + spat + log(offset);

  // ===========================================================================
  // Observation likelihood
  // ---------------------------------------------------------------------------
  Type disp = tweedie_phi(obs_lik_pars(0));
  Type shape = tweedie_p(obs_lik_pars(1));
  REPORT(disp);
  REPORT(shape);

  for (int i = 0; i < N_obs; i++) {
    jnll(1) -= dtweedie(obs(i), tweedie_mu(linpred(i)), disp, shape, true);
  }

  SIMULATE {
    for (int i = 0; i < N_obs; i++) {
      obs(i) = rtweedie(exp(linpred(i)), disp, shape);
      REPORT(obs);
    }
  }

  // ===========================================================================
  // Reports
  // ---------------------------------------------------------------------------
  vector<Type> rho_sp;
  vector<Type> sigma_sp;
  rho_sp = sqrt(8) / kappa;
  sigma_sp = tau / (kappa * 2 * sqrt(PI));

  ADREPORT(rho_sp);
  ADREPORT(sigma_sp);
  ADREPORT(beta_scale);

  return jnll.sum();
}
