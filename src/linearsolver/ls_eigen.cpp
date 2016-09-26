#ifndef _LS_EIGEN_CPP
#define _LS_EIGEN_CPP
#include "ls_eigen.h"
/*!
  \brief High level linear solver wrapper to the Eigen linear algebra library. 
*/
template <class Tx, class Tad>
LinearSolverEigen<Tx, Tad>::LinearSolverEigen(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
	const auto nic = val_mesh->nic;
	const auto njc = val_mesh->njc;
	const auto nq = val_mesh->solution->nq;
	const auto ntrans = val_mesh->solution->ntrans;
	n = nic*njc*(nq+ntrans);
		
	rhs = Eigen::MatrixXd(n, 1);
	dq = Eigen::MatrixXd(n, 1);
	lhs = Eigen::SparseMatrix<Tx, Eigen::ColMajor>(n, n);

	solver = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<Tx>>>();
};

template <class Tx, class Tad>
void LinearSolverEigen<Tx, Tad>::preallocate(int nnz){
	lhs.reserve(Eigen::VectorXi::Constant(n,MAX_NNZ));
};

template <class Tx, class Tad>
LinearSolverEigen<Tx, Tad>::~LinearSolverEigen(){
};

template <class Tx, class Tad>
void LinearSolverEigen<Tx, Tad>::set_lhs(int nnz, unsigned int *rind, unsigned int *cind, Tx *values){
	number_lhs_update += 1;
	Tx value_tmp;
	for(size_t i=0; i<nnz; i++){
		value_tmp = values[i];
		lhs.coeffRef(rind[i],cind[i]) = value_tmp;
	}
};

template <class Tx, class Tad>
void LinearSolverEigen<Tx, Tad>::set_rhs(Tx *val_rhs){
	number_rhs_update += 1;
	for(size_t i=0; i<n; i++){
		rhs(i, 0) = val_rhs[i];
	}
};


template <class Tx, class Tad>
void LinearSolverEigen<Tx, Tad>::solve(bool transpose_lhs){
	if(number_rhs_update != number_lhs_update){
		spdlog::get("console")->warn("Number of rhs update does not match with the number of lhs update!");
	}
	auto lhs_t = lhs;
	if(transpose_lhs)
		auto lhs_t = lhs.transpose();

	lhs_t.makeCompressed();
	if(!pattern_analyzed){
		solver->analyzePattern(lhs_t);
		pattern_analyzed = true;
	}
	solver->factorize(lhs_t);
	dq = solver->solve(rhs);
}
	
template <class Tx, class Tad>
void LinearSolverEigen<Tx, Tad>::solve_and_update(Tx *q, Tx under_relaxation){
	solve();
	for(size_t i=0; i<n; i++){
		q[i] = q[i] + dq(i,0)*under_relaxation;
	}
};


template <class Tx, class Tad>
void LinearSolverEigen<Tx, Tad>::get_solution(Tx *solution){
	for(size_t i=0; i<n; i++){
		solution[i] = dq(i,0);
	}
}

template <class Tx, class Tad>
void LinearSolverEigen<Tx, Tad>::reset_lhs(){
	pattern_analyzed = false;
}

template class LinearSolverEigen<double,adouble>;
#endif
