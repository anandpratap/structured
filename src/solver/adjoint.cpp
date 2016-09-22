#ifndef _ADJOINT_CPP
#define _ADJOINT_CPP
#include "adjoint.h"
#include "solution.h"
#include "eulerequation.h"
#include "linearsolver.h"
template<class Tq>
Tq objective(const Array3D<const Tq>& q){
	Tq J = 0.0;
	for(size_t i=0; i<q.extent(0); i++){
		for(size_t j=0; j<q.extent(1); j++){
			J += q[i][j][1]*q[i][j][1];
		}
	}
	return J;
}

template<class Tx, class Tad>
AdjointSolver<Tx, Tad>::AdjointSolver(){};

template<class Tx, class Tad>
AdjointSolver<Tx, Tad>::~AdjointSolver(){};

template<class Tx, class Tad>
void AdjointSolver<Tx, Tad>::solve(std::shared_ptr<Mesh<Tx,Tad>> mesh){
	auto nic = mesh->nic;
	auto njc = mesh->njc;
	auto equation = mesh->equation;
	auto solution = mesh->solution;
	auto nq = solution->nq;
	auto ntrans = solution->ntrans;
	auto psi = solution->psi;
	auto linearsolver = mesh->linearsolver;
	auto dJdq = Array3D<Tx>(nic, njc, nq+ntrans);
	spdlog::get("console")->info("Calculating adjoint");
	// calc dRdq
	trace_on(2);
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			for(size_t k=0; k<nq+ntrans; k++){
				solution->a_q[i][j][k] <<= solution->q[i][j][k];
			}
		}
	}
	equation->calc_residual(solution->a_q.const_ref(), solution->a_rhs, true);
		
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			for(size_t k=0; k<nq+ntrans; k++){
				solution->a_rhs[i][j][k] >>= solution->rhs[i][j][k];
			}
		}
	}
	trace_off();
	Tx *q_ptr = solution->q.data();
	sparse_jac(2,solution->nt,solution->nt,solution->repeat,q_ptr,&solution->nnz,&solution->rind,&solution->cind,&solution->values,solution->options);

	// transpose, see that the rind and cind are swapped!!!
	linearsolver->reset_lhs();
	linearsolver->set_lhs(solution->nnz, solution->cind, solution->rind, solution->values);
	free(solution->rind); solution->rind=nullptr;
	free(solution->cind); solution->cind=nullptr;
	free(solution->values); solution->values=nullptr;

	// calc dJdq
	trace_on(3);
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			for(size_t k=0; k<nq+ntrans; k++){
				solution->a_q[i][j][k] <<= solution->q[i][j][k];
			}
		}
	}
	auto a_J = objective<Tad>(solution->a_q.const_ref());
	Tx J;
	a_J >>= J;
	trace_off();
	spdlog::get("console")->info("Objective function: {}", J);

	gradient(3, solution->nt, solution->q.data(), dJdq.data());
	linearsolver->set_rhs(dJdq.data());
	linearsolver->solve();
	linearsolver->get_solution(psi.data());
	mesh->iomanager->write_tecplot_adjoint();
	// solve for psi
};


template<class Tx, class Tad>
void AdjointSolver<Tx, Tad>::calc_sensitivity(std::shared_ptr<Mesh<Tx,Tad>> mesh){
	auto nic = mesh->nic;
	auto njc = mesh->njc;
	auto equation = mesh->equation;
	auto solution = mesh->solution;
	auto design_parameters = mesh->design_parameters;
	auto nq = solution->nq;
	auto ntrans = solution->ntrans;
	auto psi = solution->psi;
	auto linearsolver = mesh->linearsolver;
	auto dJdbeta = Array2D<Tx>(nic, njc);
	dJdbeta.fill(0.0);

	trace_on(4);
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			design_parameters->a_beta[i][j] <<= design_parameters->beta[i][j];
		}
	}
	equation->calc_residual(solution->a_q.const_ref(), solution->a_rhs, true);
		
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			for(size_t k=0; k<nq+ntrans; k++){
				solution->a_rhs[i][j][k] >>= solution->rhs[i][j][k];
			}
		}
	}
	trace_off();
	auto q_ptr = solution->q.data();
	sparse_jac(4,solution->nt,design_parameters->n,solution->repeat,design_parameters->beta.data(),&solution->nnz,&solution->rind,&solution->cind,&solution->values,solution->options);



	auto psi_ptr = solution->psi.data();
	auto dJdbeta_ptr = dJdbeta.data();
	for(size_t i=0; i<solution->nnz; i++){
		dJdbeta_ptr[solution->cind[i]] -= psi_ptr[solution->rind[i]]*solution->values[i];
	}

	auto sum_ = 0.0;
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			sum_ += dJdbeta[i][j];
		}
	}
	spdlog::get("console")->info("Sum sens: {}", sum_);
	free(solution->rind); solution->rind=nullptr;
	free(solution->cind); solution->cind=nullptr;
	free(solution->values); solution->values=nullptr;

}

#if defined(ENABLE_ADOLC)
template class AdjointSolver<double, adouble>;
#else
#endif
#endif
