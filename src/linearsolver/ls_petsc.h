#ifndef _PETSC_H
#define _PETSC_H
#define CONFIG_PETSC_TOL 1e-12
#define CONFIG_PETSC_MAXITER 1000
#include "petsc.h"
#include "petscksp.h"

/*!
  \brief High level linear solver wrapper to the Petsc linear algebra library. 
*/
template <class Tx, class Tad>
class LinearSolverPetsc{
private:
	size_t n; //!< Size of the linear system, square matrix is assumed.
	Vec rhs; //!< Vector container for the right hand side
	Vec dq; //!< Vector container for the solution
	Mat lhs; //!< Matrix container for the left hand side
	Tx *dq_array; //!< Pointer to data of the solution
	PC pc; //!< Petsc preconditioner
	KSP ksp; //!< Petsc Krylov space solver

	size_t number_lhs_update = 0; //!< Number of times the left hand side is updated.
	size_t number_rhs_update = 0; //!< Number of times the right hand side is updated.
public:
	LinearSolverPetsc(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config);
	~LinearSolverPetsc();
	void preallocate(int nnz);	
	void set_lhs(int nnz, unsigned int *rind, unsigned int *cind, Tx *values);
	void set_rhs(Tx *val_rhs);	
	void solve();
	void solve_and_update(Tx *q, Tx under_relaxation);
};

#endif
