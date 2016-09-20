#ifndef _PETSC_H
#define _PETSC_H
#define CONFIG_PETSC_TOL 1e-12
#define CONFIG_PETSC_MAXITER 1000
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
	LinearSolverPetsc(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
		const auto nic = val_mesh->nic;
		const auto njc = val_mesh->njc;
		const auto nq = val_mesh->solution->nq;
		const auto ntrans = val_mesh->solution->ntrans;
		n = nic*njc*(nq+ntrans);
		PetscInitialize(&val_config->argc, &val_config->argv, NULL, NULL);

		PetscErrorCode ierr;
		ierr = VecCreate(PETSC_COMM_WORLD, &dq);
		ierr = PetscObjectSetName((PetscObject) dq, "Solution");
		ierr = VecSetSizes(dq, PETSC_DECIDE, n);
		ierr = VecSetFromOptions(dq);
		ierr = VecDuplicate(dq,&rhs);

		MatCreateSeqAIJ(PETSC_COMM_WORLD, n, n, MAX_NNZ, NULL, &lhs);
		
		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
		ierr = KSPSetOperators(ksp, lhs, lhs);
		ierr = KSPGetPC(ksp, &pc);
		ierr = PCSetType(pc, PCLU);
		ierr = KSPSetType(ksp, KSPGMRES);
		ierr = KSPSetTolerances(ksp, 1e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
		ierr = KSPSetFromOptions(ksp);
		MatSetBlockSize(lhs, nq+ntrans);
	};

	~LinearSolverPetsc(){
		VecDestroy(&rhs);
		VecDestroy(&dq);
		MatDestroy(&lhs);
		PetscFinalize();
	};

	void preallocate(int nnz){

	};
	
	void set_lhs(int nnz, unsigned int *rind, unsigned int *cind, Tx *values){
		number_lhs_update += 1;
		PetscScalar value_tmp;
		PetscErrorCode ierr;
		int row_idx, col_idx;
		for(size_t i=0; i<nnz; i++){
			value_tmp = values[i];
			row_idx = rind[i];
			col_idx = cind[i];
			MatSetValue(lhs, row_idx, col_idx, value_tmp ,INSERT_VALUES);
		}
		MatAssemblyBegin(lhs, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(lhs, MAT_FINAL_ASSEMBLY);
	};

	void set_rhs(Tx *val_rhs){
		number_rhs_update += 1;
		PetscScalar value_tmp;
		PetscErrorCode ierr;
		for(size_t i=0; i<n; i++){
			value_tmp = val_rhs[i];
			VecSetValue(rhs, i, value_tmp, INSERT_VALUES);
		}
	};
	
	void solve(){
		if(number_rhs_update != number_lhs_update){
			spdlog::get("console")->warn("Number of rhs update does not match with the number of lhs update!");
		}
		KSPSolve(ksp, rhs, dq);
	};

	void solve_and_update(Tx *q, Tx under_relaxation){
		solve();
		VecGetArray(dq, &dq_array);
		for(size_t i=0; i<n; i++){
			q[i] = q[i] + dq_array[i]*under_relaxation;
		}
	};   
};

#endif
