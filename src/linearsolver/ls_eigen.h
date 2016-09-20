#ifndef _EIGEN_H
#define _EIGEN_H
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
	LinearSolverEigen(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
		const auto nic = val_mesh->nic;
		const auto njc = val_mesh->njc;
		const auto nq = val_mesh->solution->nq;
		const auto ntrans = val_mesh->solution->ntrans;
		n = nic*njc*(nq+ntrans);
		
		rhs = Eigen::MatrixXd(n, 1);
		dq = Eigen::MatrixXd(n, 1);
		lhs = Eigen::SparseMatrix<Tx, Eigen::ColMajor>(n, n);

		solver = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<Tx>>>();
	};

	void preallocate(int nnz){
		lhs.reserve(Eigen::VectorXi::Constant(n,MAX_NNZ));
	};

	~LinearSolverEigen(){
	};

	void set_lhs(int nnz, unsigned int *rind, unsigned int *cind, Tx *values){
		number_lhs_update += 1;
		Tx value_tmp;
		for(size_t i=0; i<nnz; i++){
			value_tmp = values[i];
			lhs.coeffRef(rind[i],cind[i]) = value_tmp;
		}
	};

	void set_rhs(Tx *val_rhs){
		number_rhs_update += 1;
		for(size_t i=0; i<n; i++){
			rhs(i, 0) = val_rhs[i];
		}
	};


	void solve(){
		if(number_rhs_update != number_lhs_update){
			spdlog::get("console")->warn("Number of rhs update does not match with the number of lhs update!");
		}
		
		lhs.makeCompressed();
		if(!pattern_analyzed){
			solver->analyzePattern(lhs);
			pattern_analyzed = true;
		}
		solver->factorize(lhs);
		dq = solver->solve(rhs);
	}
	
	void solve_and_update(Tx *q, Tx under_relaxation){
		solve();
		for(size_t i=0; i<n; i++){
			q[i] = q[i] + dq(i,0)*under_relaxation;
		}
	};
};

#endif
