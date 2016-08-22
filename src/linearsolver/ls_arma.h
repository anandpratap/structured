#ifndef _ARMA_H
#define _ARMA_H

class LinearSolverArma{
 public:
	Mesh<double> *mesh;
	arma::mat rhs, dq;
	arma::sp_mat jac;
	LinearSolverArma(Mesh<double> *val_mesh, Config *val_config){
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
		for(uint i=0; i<nnz; i++){
			value_tmp = values[i];
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

	void solve_and_update(double *q, double UNDER_RELAXATION){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;

		dq = arma::spsolve(jac, rhs, "superlu");
		for(int i=0; i<nic*njc*nq; i++){
			q[i] = q[i] + dq(i,0)*UNDER_RELAXATION;
		}
	};

	void preallocate(int nnz){

	};
};

#endif
