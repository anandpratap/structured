#ifndef _SOLVER_CPP
#define _SOLVER_CPP
#include "solver.h"
template<class Tx>
void set_rarray(size_t size, Tx* __restrict__ dest, Tx* __restrict__ src){
#pragma omp parallel for
	for(size_t i=0; i<size; i++){
		dest[i] = src[i];
	}
}

template<class Tx>
void update_forward_euler(size_t size, Tx* __restrict__ q, Tx* __restrict__ rhs, Tx* __restrict__ dt){
#pragma omp parallel for
	for(size_t i=0; i<size; i++){
		q[i] = q[i] + rhs[i]*dt[i];
	}
}

template<class Tx, class To>
void update_rk4(size_t size, Tx* __restrict__ q_i, Tx* __restrict__ q, Tx* __restrict__ rhs, Tx* __restrict__ dt, To order){
#pragma omp parallel for
	for(size_t i=0; i<size; i++){
		q_i[i] = q[i] + rhs[i]*dt[i]/(4.0 - order);
	}
}

template<class Tx, class Tad>
void Solver<Tx, Tad>::add_mesh(std::shared_ptr<Mesh<Tx,Tad>> mesh){
		mesh_list.push_back(mesh);
}

template <class Tx, class Tad>
Solver<Tx, Tad>::Solver(std::shared_ptr<Config<Tx>> val_config){
	config = val_config;

	label = config->io->label;
	logger_convergence = spdlog::basic_logger_mt("convergence", label+".history", true);
	logger_convergence->info(" ");
	logger = spdlog::get("console");

	CFL = config->solver->cfl;
	UNDER_RELAXATION = config->solver->under_relaxation;
	

}
template <class Tx, class Tad>
Solver<Tx, Tad>::~Solver(){
}


template <class Tx, class Tad>
bool Solver<Tx, Tad>::step(std::shared_ptr<Mesh<Tx,Tad>> mesh, size_t counter, Tx t){
	auto ni = mesh->ni;
	auto nj = mesh->nj;
	auto nic = mesh->nic;
	auto njc = mesh->njc;
	auto nq = mesh->solution->nq;
	auto ntrans = mesh->solution->ntrans;
	Tx l2norm[nq+ntrans] = {1e10};

	auto iomanager = mesh->iomanager;
	auto equation = mesh->equation;
	auto solution = mesh->solution;
	
	equation->calc_dt(CFL);
	config->profiler->reset_time_residual();
	
#if defined(ENABLE_ADOLC)

		
	trace_on(1);
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

	if(config->solver->order != config->solver->lhs_order){
		equation->calc_residual(solution->a_q.const_ref(), solution->a_rhs, false);
		for(size_t i=0; i<nic; i++){
			for(size_t j=0; j<njc; j++){
				for(size_t k=0; k<nq+ntrans; k++){
					solution->rhs[i][j][k] = value(solution->a_rhs[i][j][k]);
				}
			}
		}
	}
#else
	if(config->solver->scheme == "forward_euler"){
		equation->calc_residual(solution->q.const_ref(), solution->rhs);
		update_forward_euler(solution->q.size(), solution->q.data(), solution->rhs.data(), solution->dt.data());
	}
	else if(config->solver->scheme == "rk4_jameson"){

		for(size_t order=0; order<4; order++){
			equation->calc_residual(solution->q_tmp.const_ref(), solution->rhs);
			update_rk4(solution->q.size(), solution->q_tmp.data(), solution->q.data(), solution->rhs.data(), solution->dt.data(), order);
		}

		set_rarray(solution->q.size(), solution->q.data(), solution->q_tmp.data());
			
	}
		
	else{
		logger->critical("scheme not defined.");
	}
	
#endif
	config->profiler->update_time_residual();
		
	Tx l2normq;
	for(size_t k=0; k<nq+ntrans; k++){
		l2normq = 0.0;
		for(size_t i=0; i<nic; i++){
			for(size_t j=0; j<njc; j++){
				l2normq += solution->rhs[i][j][k]*solution->rhs[i][j][k];
			}
		}
		l2norm[k] = sqrt(l2normq);
	}
		
	if(counter > config->solver->iteration_max){
		logger->info("Max iteration reached!");
		iomanager->write(counter);
		auto dt_main = config->profiler->current_time();
		logger->info("Final:: Step: {:08d} Time: {:.2e} Wall Time: {:.2e} CFL: {:.2e} Density Norm: {:.2e}", counter, t, dt_main, CFL, l2norm[0]);
		logger_convergence->info("{:08d} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}", counter, t, dt_main, CFL, l2norm[0], l2norm[1], l2norm[2], l2norm[3]);
		return true;
	}
	//		logger->debug("transport l2 {}", l2norm[4]);
	if(l2norm[0] < config->solver->tolerance && l2norm[1] < config->solver->tolerance && 0){
		logger->info("Convergence reached!");
		iomanager->write(counter);
		auto dt_main = config->profiler->current_time();
		logger->info("Final:: Step: {:08d} Time: {:.2e} Wall Time: {:.2e} CFL: {:.2e} Density Norm: {:.2e}", counter, t, dt_main, CFL, l2norm[0]);
		logger_convergence->info("{:08d} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}", counter, t, dt_main, CFL, l2norm[0], l2norm[1], l2norm[2], l2norm[3]);
		return true;
	}
#if defined(ENABLE_ADOLC)
	config->profiler->reset_time_jacobian();
	Tx *q_ptr = solution->q.data();
	sparse_jac(1,solution->nt,solution->nt,solution->repeat,q_ptr,&solution->nnz,&solution->rind,&solution->cind,&solution->values,solution->options);
	config->profiler->update_time_jacobian();
	logger->debug("NNZ = {}", solution->nnz);
		
	if(counter == 0)
		mesh->linearsolver->preallocate(solution->nnz);
	Tx dt_local = 0.0;
	for(size_t i=0; i<solution->nnz; i++){
		const auto k_idx = solution->rind[i]%(nq+ntrans);
		const auto j_idx = ((solution->rind[i]-k_idx)/(nq+ntrans))%njc;
		const auto i_idx = (solution->rind[i] - k_idx - j_idx*(nq+ntrans))/njc/(nq+ntrans);
		dt_local = solution->dt[i_idx][j_idx][k_idx];
		solution->values[i] = -solution->values[i];
		//			std::cout<<dt_local<<std::endl;
		if(solution->rind[i] == solution->cind[i]){solution->values[i] += 1.0/dt_local;}
	}
	config->profiler->reset_time_linearsolver();
	mesh->linearsolver->set_lhs(solution->nnz, solution->rind, solution->cind, solution->values);
	mesh->linearsolver->set_rhs(solution->rhs.data());
	mesh->linearsolver->solve_and_update(solution->q.data(), UNDER_RELAXATION);
	auto dt_perf = config->profiler->update_time_linearsolver();
	//q[i][j][k] = q[i][j][k] + solution->rhs[i][j][k]*dt;
		
	logger->info("Linear algebra time = {:03.2f}", dt_perf);
		
	free(solution->rind); solution->rind=nullptr;
	free(solution->cind); solution->cind=nullptr;
	free(solution->values); solution->values=nullptr;
#else
#endif
	if(counter % config->io->stdout_frequency == 0){
		auto dt_main = config->profiler->current_time();
		logger->info("Step: {:08d} Time: {:.2e} Wall Time: {:.2e} CFL: {:.2e} Density Norm: {:.2e}", counter, t, dt_main, CFL, l2norm[0]);
		logger_convergence->info("{:08d} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}", counter, t, dt_main, CFL, l2norm[0], l2norm[1], l2norm[2], l2norm[3]);
	}
		
	if(counter % config->io->fileout_frequency == 0){
		iomanager->write(counter);
	}
	return false;
}
template <class Tx, class Tad>
void Solver<Tx, Tad>::solve(){
	size_t counter = 0;
	Tx t = 0.0;
	logger->info("Welcome to structured!");
	config->profiler->timer_main->reset();
	while(1){
		bool if_break = true;
		for(auto&& mesh : mesh_list)
			if_break = step(mesh, counter, t) && if_break;
		if(if_break){
			break;
		}
		t += 0;
		counter += 1;
		auto cfl_ramp =  config->solver->cfl_ramp;
		if(cfl_ramp){
			auto cfl_ramp_iteration = config->solver->cfl_ramp_iteration;
			if(counter > cfl_ramp_iteration){
				auto cfl_ramp_exponent = config->solver->cfl_ramp_exponent;
				CFL = pow(CFL, cfl_ramp_exponent);
				CFL = std::min(CFL, static_cast<Tx>(1e12));
			}
		}

		auto under_relaxation_ramp =  config->solver->under_relaxation_ramp;
		if(under_relaxation_ramp){
			auto under_relaxation_ramp_iteration = config->solver->under_relaxation_ramp_iteration;
			if(counter > under_relaxation_ramp_iteration){
				auto under_relaxation_ramp_exponent = config->solver->under_relaxation_ramp_exponent;
				UNDER_RELAXATION = pow(UNDER_RELAXATION, under_relaxation_ramp_exponent);
				UNDER_RELAXATION = std::min(UNDER_RELAXATION, static_cast<Tx>(10.0));
			}
		}


		
		
	}
}

#if defined(ENABLE_ADOLC)
template class Solver<double, adouble>;
#else
template class Solver<double, double>;
template void set_rarray<double>;
template void update_forward_euler<double>;
template void update_rk4<double,size_t>;
template class Solver<float, float>;
template void set_rarray<float>;
template void update_forward_euler<float>;
template void update_rk4<float,size_t>;
#endif
#endif
