#ifndef _ARMA_H
#define _ARMA_H

/*!
  \brief High level linear solver wrapper to the Armadillo linear algebra library. 
 */
template<class Tx, class Tad>
class LinearSolverArma{
 public:
	uint n; //!< Size of the linear system, square matrix is assumed.
	arma::Mat<Tx> rhs; //!< Matrix container for the right hand side
	arma::Mat<Tx> dq; //!< Matrix container for the solution
	arma::SpMat<Tx> lhs; //!< Matrix container for the left hand side

	uint number_lhs_update = 0; //!< Number of times the left hand side is updated.
	uint number_rhs_update = 0; //!< Number of times the right hand side is updated.

	LinearSolverArma(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
		const auto nic = val_mesh->nic;
		const auto njc = val_mesh->njc;
		const auto nq = val_mesh->solution->nq;
		const auto ntrans = val_mesh->solution->ntrans;
		n = nic*njc*(nq+ntrans);
		rhs = arma::Mat<Tx>(n, 1);
		dq = arma::Mat<Tx>(n, 1);
		lhs = arma::SpMat<Tx>(n, n);
	};

	void set_lhs(int nnz, unsigned int *rind, unsigned int *cind, Tx *values){
		number_lhs_update += 1;
		Tx value_tmp;
		for(uint i=0; i<nnz; i++){
			value_tmp = values[i];
			lhs(rind[i],cind[i]) = value_tmp;
		}
	};

	void set_rhs(Tx *val_rhs){
		number_rhs_update += 1;
		for(uint i=0; i<n; i++){
			rhs(i, 0) = val_rhs[i];
		}
	};

	void solve(){
		if(number_rhs_update != number_lhs_update){
			spdlog::get("console")->warn("Number of rhs update does not match with the number of lhs update!");
		}
		dq = arma::spsolve(lhs, rhs, "superlu");
	};
	
	void solve_and_update(Tx *q, Tx under_relaxation){
		for(int i=0; i<n; i++){
			q[i] = q[i] + dq(i,0)*under_relaxation;
		}
	};

	void preallocate(int nnz){

	};
};

#endif
