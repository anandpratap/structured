#ifndef _LINEARSOLVER_H
#define _LINEARSOLVER_H
#include "common.h"


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
		double dt = 1e2;
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

		dq = arma::spsolve(jac, rhs, "lapack");
		for(int i=0; i<nic*njc*nq; i++){
			q[i] = q[i] + dq(i,0);
		}
	};

	void preallocate(int nnz){

	};
};

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
		double dt = 1e2;
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
			q[i] = q[i] + dq(i,0);
		}
	};

};



#endif
