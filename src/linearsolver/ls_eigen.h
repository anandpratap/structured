#ifndef _LS_EIGEN_H
#define _LS_EIGEN_H
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "mesh.h"
#include "config.h"
#include "solution.h"
#include "common.h"
/*!
  \brief High level linear solver wrapper to the Eigen linear algebra library. 
*/
template <class Tx, class Tad>
class LinearSolverEigen{
private:
	size_t n; //!< Size of the linear system, square matrix is assumed.
	Eigen::MatrixXd rhs; //!< Matrix container for the right hand side
	Eigen::MatrixXd dq; //!< Matrix container for the solution
	Eigen::SparseMatrix<Tx,  Eigen::ColMajor> lhs; //!< Matrix container for the left hand side
	std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<Tx>>> solver; //!< Pointer to the linear solver

	bool pattern_analyzed = false; //!< To ensure that pattern is analyzed when the solver is called for the first time 
	size_t number_lhs_update = 0; //!< Number of times the left hand side is updated.
	size_t number_rhs_update = 0; //!< Number of times the right hand side is updated.
public:
	LinearSolverEigen(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config);
	void preallocate(int nnz);
	~LinearSolverEigen();
	void set_lhs(int nnz, unsigned int *rind, unsigned int *cind, Tx *values);
	void set_rhs(Tx *val_rhs);
	void solve();	
	void solve_and_update(Tx *q, Tx under_relaxation);
	void get_solution(Tx *solution);
	void reset_lhs();
};
#endif
