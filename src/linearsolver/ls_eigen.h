#ifndef _EIGEN_H
#define _EIGEN_H

class LinearSolverEigen{
	Mesh<double> *mesh;
	Eigen::MatrixXd rhs, dq;
	Eigen::SparseMatrix<double, Eigen::ColMajor> jac;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> *solver;
	unsigned int c = 0;
 public:
	LinearSolverEigen(Mesh<double> *val_mesh, Config *val_config){
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
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;
		jac.reserve(Eigen::VectorXi::Constant(nic*njc*nq,36));
	};

	~LinearSolverEigen(){
		delete solver;
	};

	void set_jac(int nnz, unsigned int *rind, unsigned int *cind, double *values){
		double value_tmp;
		for(uint i=0; i<nnz; i++){
			value_tmp = values[i];
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

	void solve_and_update(double *q, double UNDER_RELAXATION){
		uint nic = mesh->nic;
		uint njc = mesh->njc;
		uint nq = mesh->solution->nq;
		jac.makeCompressed();
		if(c == 0){
			solver->analyzePattern(jac);
			c = 1;
		}
		solver->factorize(jac);
		dq = solver->solve(rhs);
		for(int i=0; i<nic*njc*nq; i++){
			q[i] = q[i] + dq(i,0)*UNDER_RELAXATION;
		}
	};

};

#endif
