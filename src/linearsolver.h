#ifndef _LINEARSOLVER_H
#define _LINEARSOLVER_H
#include "common.h"
double dt_global = 1e1;
double fac_global = 0.6;
#if defined(ENABLE_ARMA)
class LinearSolverArma{
 public:
	Mesh<double> *mesh;
	arma::mat rhs, dq;
	arma::sp_mat jac;
	LinearSolverArma(Mesh<double> *val_mesh){
		mesh = val_mesh;

		uint ni = mesh->ni;
		uint nj = mesh->nj;
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;

		rhs = arma::mat(nic*njc*nq, 1);
		dq = arma::mat(nic*njc*nq, 1);
		jac = arma::sp_mat(nic*njc*nq, nic*njc*nq);
	};

	void set_jac(int nnz, unsigned int *rind, unsigned int *cind, double *values){
		double value_tmp;
		double dt = dt_global;
		for(uint i=0; i<nnz; i++){
			value_tmp = -values[i];
			if(rind[i] == cind[i]){value_tmp += 1.0/dt;}
			jac(rind[i],cind[i]) = value_tmp;
		}
	};

	void set_rhs(double *val_rhs){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;

		for(uint i=0; i<nic*njc*nq; i++){
			rhs(i, 0) = val_rhs[i];
		}
	};

	void solve_and_update(double *q){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;

		dq = arma::spsolve(jac, rhs, "superlu");
		for(int i=0; i<nic*njc*nq; i++){
			q[i] = q[i] + dq(i,0)*fac_global;
		}
	};

	void preallocate(int nnz){

	};
};

#endif

#if defined(ENABLE_EIGEN)
class LinearSolverEigen{
	Mesh<double> *mesh;
	Eigen::MatrixXd rhs, dq;
	Eigen::SparseMatrix<double, Eigen::ColMajor> jac;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> *solver;

 public:
	LinearSolverEigen(Mesh<double> *val_mesh){
		mesh = val_mesh;

		uint ni = mesh->ni;
		uint nj = mesh->nj;
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;
	
		rhs = Eigen::MatrixXd(nic*njc*nq, 1);
		dq = Eigen::MatrixXd(nic*njc*nq, 1);
		jac = Eigen::SparseMatrix<double,Eigen::ColMajor>(nic*njc*nq, nic*njc*nq);

		solver = new Eigen::SparseLU<Eigen::SparseMatrix<double>>();
	};

	void preallocate(int nnz){
		jac.reserve(nnz);
	};

	~LinearSolverEigen(){
		delete solver;
	};

	void set_jac(int nnz, unsigned int *rind, unsigned int *cind, double *values){
		double value_tmp;
		double dt = dt_global;
		for(uint i=0; i<nnz; i++){
			value_tmp = -values[i];
			if(rind[i] == cind[i]){value_tmp += 1.0/dt;}
			jac.coeffRef(rind[i],cind[i]) = value_tmp;
		}
	};

	void set_rhs(double *val_rhs){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;

		for(uint i=0; i<nic*njc*nq; i++){
			rhs(i, 0) = val_rhs[i];
		}
	};

	void solve_and_update(double *q){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;
		jac.makeCompressed();
		solver->analyzePattern(jac);
		solver->factorize(jac);
		dq = solver->solve(rhs);
		for(int i=0; i<nic*njc*nq; i++){
			q[i] = q[i] + dq(i,0)*fac_global;
		}
	};

};

#endif

#if defined(ENABLE_PETSC)
#define CONFIG_PETSC_TOL 1e-12
#define CONFIG_PETSC_MAXITER 1000

class LinearSolverPetsc{
	Mesh<double> *mesh;
	Vec dq, rhs;
	Mat jac;
	double *dq_array;
	PC pc;
	KSP ksp;

 public:
	LinearSolverPetsc(Mesh<double> *val_mesh){
		mesh = val_mesh;
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;
		int nvar = nic*njc*nq;
		PetscInitialize(NULL, NULL, NULL, NULL);

		PetscErrorCode ierr;
		ierr = VecCreate(PETSC_COMM_WORLD, &dq);
		ierr = PetscObjectSetName((PetscObject) dq, "Solution");
		ierr = VecSetSizes(dq, PETSC_DECIDE, nvar);
		ierr = VecSetFromOptions(dq);
		ierr = VecDuplicate(dq,&rhs);

		MatCreateSeqAIJ(PETSC_COMM_WORLD, nvar, nvar, 36, NULL, &jac);

		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
		ierr = KSPSetOperators(ksp, jac, jac);
		ierr = KSPGetPC(ksp, &pc);
		ierr = PCSetType(pc, PCLU);
		ierr = KSPSetTolerances(ksp, CONFIG_PETSC_TOL, PETSC_DEFAULT, PETSC_DEFAULT, CONFIG_PETSC_MAXITER);
		ierr = KSPSetFromOptions(ksp);
		KSPSetType(ksp, KSPGMRES);
			
	};

	~LinearSolverPetsc(){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;
		VecDestroy(&rhs);
		VecDestroy(&dq);
		MatDestroy(&jac);
		PetscFinalize();
	};

	void preallocate(int nnz){

	};
	
	void set_jac(int nnz, unsigned int *rind, unsigned int *cind, double *values){
		PetscScalar value_tmp;
		PetscErrorCode ierr;
		int row_idx, col_idx;
		double dt = dt_global;
		for(uint i=0; i<nnz; i++){
			value_tmp = -values[i];
			if(rind[i] == cind[i]){value_tmp += 1.0/dt;}
			row_idx = rind[i];
			col_idx = cind[i];
			MatSetValue(jac, row_idx, col_idx, value_tmp ,INSERT_VALUES);
		}
		MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
	};

	void set_rhs(double *val_rhs){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;
		PetscScalar value_tmp;
		PetscErrorCode ierr;

		for(int i=0; i<nic*njc*nq; i++){
			value_tmp = val_rhs[i];
			VecSetValue(rhs, i, value_tmp, INSERT_VALUES);
		}
	};

	void solve_and_update(double *q){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;
		KSPSolve(ksp, rhs, dq);
		VecGetArray(dq, &dq_array);

		for(int i=0; i<nic*njc*nq; i++){
			q[i] = q[i] + dq_array[i]*fac_global;
		}
	};   
};

#endif

#endif
