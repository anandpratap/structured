#ifndef _LS_ARMA_H
#define _LS_ARMA_H
#include <armadillo>
#include "common.h"
#include "mesh.h"
#include "config.h"
#include "solution.h"

/*!
  \brief High level linear solver wrapper to the Armadillo linear algebra library. 
*/
template<class Tx, class Tad>
class LinearSolverArma{
public:
	size_t n; //!< Size of the linear system, square matrix is assumed.
	arma::Mat<Tx> rhs; //!< Matrix container for the right hand side
	arma::Mat<Tx> dq; //!< Matrix container for the solution
	arma::SpMat<Tx> lhs; //!< Matrix container for the left hand side

	size_t number_lhs_update = 0; //!< Number of times the left hand side is updated.
	size_t number_rhs_update = 0; //!< Number of times the right hand side is updated.

	LinearSolverArma(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config);
	void set_lhs(int nnz, unsigned int *rind, unsigned int *cind, Tx *values);
	void set_rhs(Tx *val_rhs);
	void solve();	
	void solve_and_update(Tx *q, Tx under_relaxation);
	void preallocate(int nnz);
};

#endif
