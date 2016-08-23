#ifndef _EIGEN_H
#define _EIGEN_H

template <class T>
class LinearSolverEigen{
	std::shared_ptr<Mesh<T>> mesh;
	Eigen::MatrixXd rhs, dq;
	Eigen::SparseMatrix<T, Eigen::ColMajor> jac;
	std::shared_ptr<Eigen::SparseLU<Eigen::SparseMatrix<T>>> solver;
	uint c = 0;
 public:
	LinearSolverEigen(std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config<T>> val_config){
		mesh = val_mesh;

		const auto ni = mesh->ni;
		const auto nj = mesh->nj;
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
	
		rhs = Eigen::MatrixXd(nic*njc*nq, 1);
		dq = Eigen::MatrixXd(nic*njc*nq, 1);
		jac = Eigen::SparseMatrix<T,Eigen::ColMajor>(nic*njc*nq, nic*njc*nq);

		solver = std::make_shared<Eigen::SparseLU<Eigen::SparseMatrix<T>>>();
	};

	void preallocate(int nnz){
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		jac.reserve(Eigen::VectorXi::Constant(nic*njc*nq,36));
	};

	~LinearSolverEigen(){
	};

	void set_jac(int nnz, unsigned int *rind, unsigned int *cind, T *values){
		T value_tmp;
		for(uint i=0; i<nnz; i++){
			value_tmp = values[i];
			jac.coeffRef(rind[i],cind[i]) = value_tmp;
		}
	};

	void set_rhs(T *val_rhs){
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;

		for(uint i=0; i<nic*njc*nq; i++){
			rhs(i, 0) = val_rhs[i];
		}
	};

	void solve_and_update(T *q, T under_relaxation){
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		jac.makeCompressed();
		if(c == 0){
			solver->analyzePattern(jac);
			c = 1;
		}
		solver->factorize(jac);
		dq = solver->solve(rhs);
		for(int i=0; i<nic*njc*nq; i++){
			q[i] = q[i] + dq(i,0)*under_relaxation;
		}
	};

};

#endif
