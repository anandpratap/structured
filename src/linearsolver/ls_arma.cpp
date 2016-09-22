#ifndef _LS_ARMA_CPP
#define _LS_ARMA_CPP
#include "ls_arma.h"
/*!
  \brief High level linear solver wrapper to the Armadillo linear algebra library. 
*/
template<class Tx, class Tad>
LinearSolverArma<Tx, Tad>::LinearSolverArma(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
		const auto nic = val_mesh->nic;
		const auto njc = val_mesh->njc;
		const auto nq = val_mesh->solution->nq;
		const auto ntrans = val_mesh->solution->ntrans;
		n = nic*njc*(nq+ntrans);
		rhs = arma::Mat<Tx>(n, 1);
		dq = arma::Mat<Tx>(n, 1);
		lhs = arma::SpMat<Tx>(n, n);
};

template<class Tx, class Tad>
void LinearSolverArma<Tx, Tad>::set_lhs(int nnz, unsigned int *rind, unsigned int *cind, Tx *values){
	number_lhs_update += 1;
		Tx value_tmp;
		for(size_t i=0; i<nnz; i++){
			value_tmp = values[i];
			lhs(rind[i],cind[i]) = value_tmp;
		}
};
template<class Tx, class Tad>
void LinearSolverArma<Tx, Tad>::set_rhs(Tx *val_rhs){
	number_rhs_update += 1;
	for(size_t i=0; i<n; i++){
		rhs(i, 0) = val_rhs[i];
	}
};

template<class Tx, class Tad>
void LinearSolverArma<Tx, Tad>::solve(){
		if(number_rhs_update != number_lhs_update){
			spdlog::get("console")->warn("Number of rhs update does not match with the number of lhs update!");
		}
		dq = arma::spsolve(lhs, rhs, "superlu");
	};
	
template<class Tx, class Tad>
void LinearSolverArma<Tx, Tad>::solve_and_update(Tx *q, Tx under_relaxation){
	for(size_t i=0; i<n; i++){
		q[i] = q[i] + dq(i,0)*under_relaxation;
	}
};

template<class Tx, class Tad>
void LinearSolverArma<Tx, Tad>::preallocate(int nnz){

};


template class LinearSolverArma<double,adouble>;
#endif
