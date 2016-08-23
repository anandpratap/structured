#ifndef _ARMA_H
#define _ARMA_H
template<class T>
class LinearSolverArma{
 public:
	std::shared_ptr<Mesh<T>> mesh;
	arma::Mat<T> rhs, dq;
	arma::SpMat<T> jac;
	LinearSolverArma(std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config> val_config){
		mesh = val_mesh;

		uint ni = mesh->ni;
		uint nj = mesh->nj;
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;

		rhs = arma::Mat<T>(nic*njc*nq, 1);
		dq = arma::Mat<T>(nic*njc*nq, 1);
		jac = arma::SpMat<T>(nic*njc*nq, nic*njc*nq);
	};

	void set_jac(int nnz, unsigned int *rind, unsigned int *cind, T *values){
		T value_tmp;
		for(uint i=0; i<nnz; i++){
			value_tmp = values[i];
			jac(rind[i],cind[i]) = value_tmp;
		}
	};

	void set_rhs(T *val_rhs){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;

		for(uint i=0; i<nic*njc*nq; i++){
			rhs(i, 0) = val_rhs[i];
		}
	};

	void solve_and_update(T *q, T under_relaxation){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;

		dq = arma::spsolve(jac, rhs, "superlu");
		for(int i=0; i<nic*njc*nq; i++){
			q[i] = q[i] + dq(i,0)*under_relaxation;
		}
	};

	void preallocate(int nnz){

	};
};

#endif
