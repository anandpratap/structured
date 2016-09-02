#ifndef _ARMA_H
#define _ARMA_H
template<class Tx>
class LinearSolverArma{
 public:
	std::shared_ptr<Mesh<Tx>> mesh;
	arma::Mat<Tx> rhs, dq;
	arma::SpMat<Tx> jac;
	LinearSolverArma(std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
		mesh = val_mesh;

		const auto ni = mesh->ni;
		const auto nj = mesh->nj;
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ntrans = mesh->solution->ntrans;
		rhs = arma::Mat<Tx>(nic*njc*(nq+ntrans), 1);
		dq = arma::Mat<Tx>(nic*njc*(nq+ntrans), 1);
		jac = arma::SpMat<Tx>(nic*njc*(nq+ntrans), nic*njc*(nq+ntrans));
	};

	void set_jac(int nnz, unsigned int *rind, unsigned int *cind, Tx *values){
		Tx value_tmp;
		for(uint i=0; i<nnz; i++){
			value_tmp = values[i];
			jac(rind[i],cind[i]) = value_tmp;
		}
	};

	void set_rhs(Tx *val_rhs){
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ntrans = mesh->solution->ntrans;
		for(uint i=0; i<nic*njc*(nq+ntrans); i++){
			rhs(i, 0) = val_rhs[i];
		}
	};

	void solve_and_update(Tx *q, Tx under_relaxation){
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ntrans = mesh->solution->ntrans;

		dq = arma::spsolve(jac, rhs, "superlu");
		for(int i=0; i<nic*njc*(nq+ntrans); i++){
			q[i] = q[i] + dq(i,0)*under_relaxation;
		}
	};

	void preallocate(int nnz){

	};
};

#endif
