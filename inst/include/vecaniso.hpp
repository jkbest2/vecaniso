using namespace Eigen;
using namespace tmbutils;
using namespace density;
using namespace R_inla;

/* Construct H anisotropy matrix in place */
template<class Type>
matrix<Type> aniso_h(matrix<Type> &H, Type theta, Type beta) {
    vector<Type> v(2);
    v(0) = cos(theta);
    v(1) = sin(theta);

    Type gamma = Type(0.5) * (sqrt(beta * beta + Type(4)) - beta);

    // matrix<Type> H(2, 2);
    H(0, 0) = gamma + v(0) * v(0);
    H(0, 1) = v(0) * v(1);
    H(1, 0) = H(0, 1);
    H(1, 1) = gamma + v(1) * v(1);

    return H;
}

/* Construct adjugate of H anisotropy matrix in place */
template<class Type>
matrix<Type> aniso_adj_h(matrix<Type> &adj_H, Type theta, Type beta) {
    vector<Type> v(2);
    v(0) = cos(theta);
    v(1) = sin(theta);

    Type gamma = Type(0.5) * (sqrt(beta * beta + Type(4)) - beta);

    adj_H(0, 0) = gamma + beta * v(1) * v(1);
    adj_H(0, 1) = -beta * v(0) * v(1);
    adj_H(1, 0) = adj_H(0, 1);
    adj_H(1, 1) = gamma + beta * v(0) * v(0);

    return adj_H;
}


/** \brief Object containing all elements of an vector-anisotropic SPDE object */
template<class Type>
struct spde_vecaniso_t{
    int n_s;
    int n_tri;
    vector<Type> Tri_Area;
    matrix<Type> E0;
    matrix<Type> E1;
    matrix<Type> E2;
    matrix<int>  TV;
    SparseMatrix<Type> G0;
    SparseMatrix<Type> G0_inv;
    vector<Type> theta;
    vector<Type> beta;

    spde_vecaniso_t(SEXP x){  /* x = List passed from R */
        n_s = 	CppAD::Integer(asVector<Type>(getListElement(x,"n_s"))[0]);
        n_tri = 	CppAD::Integer(asVector<Type>(getListElement(x,"n_tri"))[0]);
        Tri_Area = asVector<Type>(getListElement(x,"Tri_Area"));
        E0 = asMatrix<Type>(getListElement(x,"E0"));
        E1 = asMatrix<Type>(getListElement(x,"E1"));
        E2 = asMatrix<Type>(getListElement(x,"E2"));
        TV = asMatrix<int>(getListElement(x,"TV"));
        G0 = asSparseMatrix<Type>(getListElement(x,"G0"));
        G0_inv = asSparseMatrix<Type>(getListElement(x,"G0_inv"));
        theta = asVector<Type>(getListElement(x, "theta"));
        beta = asVector<Type>(getListElement(x, "beta"));
    }
};


/** Modified from TMB/inst/include/tmbutils/R_inla.hpp **/
template<class Type>
SparseMatrix<Type> Q_spde(spde_vecaniso_t<Type> spde, Type kappa, Type beta_scale){

    int i;
    Type kappa_pow2 = kappa*kappa;
    Type kappa_pow4 = kappa_pow2*kappa_pow2;

    int n_s = spde.n_s;
    int n_tri = spde.n_tri;
    vector<Type> Tri_Area = spde.Tri_Area;
    matrix<Type> E0 = spde.E0;
    matrix<Type> E1 = spde.E1;
    matrix<Type> E2 = spde.E2;
    matrix<int> TV = spde.TV;
    SparseMatrix<Type> G0 = spde.G0;
    SparseMatrix<Type> G0_inv = spde.G0_inv;
    vector<Type> theta = spde.theta;
    vector<Type> beta = spde.beta;

    //Type H_trace = H(0,0)+H(1,1);
    //Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);
    SparseMatrix<Type> G1_aniso(n_s,n_s);
    SparseMatrix<Type> G2_aniso(n_s,n_s);
    // Calculate adjugate of H
    matrix<Type> adj_H(2,2);
    // adj_H(0,0) = H(1,1);
    // adj_H(0,1) = -1 * H(0,1);
    // adj_H(1,0) = -1 * H(1,0);
    // adj_H(1,1) = H(0,0);
    // Calculate new SPDE matrices

    // Calculate G1 - pt. 1
    array<Type> Gtmp(n_tri,3,3);
    for(i=0; i<n_tri; i++){
        aniso_adj_h(adj_H, theta(i), beta_scale * beta(i));
        // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.
        Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
        Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
        Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
        Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
        Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
        Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    }
    // Calculate G1 - pt. 2
    for(i=0; i<n_tri; i++){
        G1_aniso.coeffRef(TV(i,1),TV(i,0)) = G1_aniso.coeffRef(TV(i,1),TV(i,0)) + (Gtmp(i,0,1));
        G1_aniso.coeffRef(TV(i,0),TV(i,1)) = G1_aniso.coeffRef(TV(i,0),TV(i,1)) + (Gtmp(i,0,1));
        G1_aniso.coeffRef(TV(i,2),TV(i,1)) = G1_aniso.coeffRef(TV(i,2),TV(i,1)) + (Gtmp(i,1,2));
        G1_aniso.coeffRef(TV(i,1),TV(i,2)) = G1_aniso.coeffRef(TV(i,1),TV(i,2)) + (Gtmp(i,1,2));
        G1_aniso.coeffRef(TV(i,2),TV(i,0)) = G1_aniso.coeffRef(TV(i,2),TV(i,0)) + (Gtmp(i,0,2));
        G1_aniso.coeffRef(TV(i,0),TV(i,2)) = G1_aniso.coeffRef(TV(i,0),TV(i,2)) + (Gtmp(i,0,2));
        G1_aniso.coeffRef(TV(i,0),TV(i,0)) = G1_aniso.coeffRef(TV(i,0),TV(i,0)) + (Gtmp(i,0,0));
        G1_aniso.coeffRef(TV(i,1),TV(i,1)) = G1_aniso.coeffRef(TV(i,1),TV(i,1)) + (Gtmp(i,1,1));
        G1_aniso.coeffRef(TV(i,2),TV(i,2)) = G1_aniso.coeffRef(TV(i,2),TV(i,2)) + (Gtmp(i,2,2));
    }
    G2_aniso = G1_aniso * G0_inv * G1_aniso;

    return kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1_aniso + G2_aniso;
}
