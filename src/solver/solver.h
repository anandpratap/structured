#ifndef _SOLVER_H
#define _SOLVER_H
#include "common.h"
#include "mesh.h"
#include "linearsolver.h"
#include "config.h"
#include "io.h"
#include "eulerequation.h"

template<class Tx>
void set_rarray(int size, Tx* __restrict__ dest, Tx* __restrict__ src){
#pragma omp parallel for
	for(auto i=0; i<size; i++){
		dest[i] = src[i];
	}
}

template<class Tx>
void update_forward_euler(int size, Tx* __restrict__ q, Tx* __restrict__ rhs, Tx* __restrict__ dt){
#pragma omp parallel for
	for(auto i=0; i<size; i++){
		q[i] = q[i] + rhs[i]*dt[i];
	}
}

template<class Tx, class To>
	void update_rk4(int size, Tx* __restrict__ q_i, Tx* __restrict__ q, Tx* __restrict__ rhs, Tx* __restrict__ dt, To order){
#pragma omp parallel for
	for(auto i=0; i<size; i++){
		q_i[i] = q[i] + rhs[i]*dt[i]/(4.0 - order);
	}
}

template<class Tx, class Tad>
class Solver{
public:
	std::shared_ptr<Mesh<Tx>> mesh;
	void solve();
	Tx UNDER_RELAXATION;
	Tx CFL;
	std::shared_ptr<Config<Tx>> config;
	std::string label;
	std::shared_ptr<spdlog::logger> logger;
	std::shared_ptr<spdlog::logger> logger_convergence;	
	std::shared_ptr<IOManager<Tx>> iomanager;

	std::shared_ptr<EulerEquation<Tx, Tad>> equation;

#if defined(ENABLE_ARMA)
	std::shared_ptr<LinearSolverArma<Tx>> linearsolver;
#endif
#if defined(ENABLE_EIGEN)
	std::shared_ptr<LinearSolverEigen<Tx>> linearsolver;
#endif
#if defined(ENABLE_PETSC)
	std::shared_ptr<LinearSolverPetsc<Tx>> linearsolver;
#endif

	
	Solver(std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> config);
	~Solver();
	void calc_dt();
	void initialize();
};
template <class Tx, class Tad>
void Solver<Tx, Tad>::calc_dt(){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;

#pragma omp parallel for
	for(uint i=0; i<mesh->nic; i++){
		for(uint j=0; j<mesh->njc; j++){
			Tx rho = mesh->solution->q[i][j][0];
			Tx u = mesh->solution->q[i][j][1]/rho;
			Tx v = mesh->solution->q[i][j][2]/rho;
			Tx rhoE = mesh->solution->q[i][j][3];
			Tx p = (rhoE - 0.5*rho*(u*u + v*v))*(GAMMA-1.0);
			Tx lambda = sqrt(GAMMA*p/rho) + fabs(u) + fabs(v);
			Tx len_min = std::min(mesh->ds_eta[i][j], mesh->ds_chi[i][j]);
			Tx mu = config->freestream->mu_inf;
			for(int k=0; k<nq; k++)
				mesh->solution->dt[i][j][k] = CFL/(lambda/len_min + 2.0*mu/len_min/len_min);
			//dt[i][j][4] = dt[i][j][0];
		}
	}

}
template <class Tx, class Tad>
void Solver<Tx, Tad>::initialize(){
	uint nic = mesh->nic;
	uint njc = mesh->njc;

	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			auto rho_inf = config->freestream->rho_inf;
			auto u_inf = config->freestream->u_inf;
			auto v_inf = config->freestream->v_inf;
			auto p_inf = config->freestream->p_inf;
			mesh->solution->q[i][j][0] = rho_inf;
			mesh->solution->q[i][j][1] = rho_inf*u_inf;
			mesh->solution->q[i][j][2] = rho_inf*v_inf;
			mesh->solution->q[i][j][3] = p_inf/(GAMMA-1.0) + 0.5*rho_inf*(u_inf*u_inf + v_inf*v_inf);
		}
	}
	if(config->io->restart){
		iomanager->read_restart();
	}
	uint nq = mesh->solution->nq;
	uint ntrans = mesh->solution->ntrans;
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			for(uint k=0; k<nq+ntrans; k++){
				mesh->solution->q_tmp[i][j][k] = mesh->solution->q[i][j][k];
			}
		}
	}
	
}

template <class Tx, class Tad>
Solver<Tx, Tad>::Solver(std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
	config = val_config;
	mesh = val_mesh;
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	uint ntrans = mesh->solution->ntrans;

#if defined(ENABLE_ARMA)
	linearsolver = std::make_shared<LinearSolverArma<Tx>>(mesh, config);
#endif

#if defined(ENABLE_EIGEN)
	linearsolver = std::make_shared<LinearSolverEigen<Tx>>(mesh, config);
#endif

#if defined(ENABLE_PETSC)
	linearsolver = std::make_shared<LinearSolverPetsc<Tx>>(mesh, config);
#endif


	label = config->io->label;
	logger_convergence = spdlog::basic_logger_mt("convergence", label+".history", true);
	logger_convergence->info(" ");
	logger = spdlog::get("console");

	CFL = config->solver->cfl;
	UNDER_RELAXATION = config->solver->under_relaxation;
	

	equation = std::make_shared<EulerEquation<Tx, Tad>>(mesh, config);
	iomanager = std::make_shared<IOManager<Tx>>(mesh, config);
}
template <class Tx, class Tad>
Solver<Tx, Tad>::~Solver(){
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	
}


template <class Tx, class Tad>
void Solver<Tx, Tad>::solve(){
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	uint ntrans = mesh->solution->ntrans;

	Tx l2norm[nq+ntrans] = {1e10};

	uint counter = 0;
	Tx t = 0.0;
	initialize();
	logger->info("Welcome to structured!");


	config->profiler->timer_main->reset();

	while(1){
		calc_dt();
		config->profiler->reset_time_residual();

#if defined(ENABLE_ADOLC)

		
		trace_on(1);
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq+ntrans; k++){
					mesh->solution->a_q[i][j][k] <<= mesh->solution->q[i][j][k];
				}
			}
		}
		equation->calc_residual(mesh->solution->a_q, mesh->solution->a_rhs);
		
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq+ntrans; k++){
					mesh->solution->a_rhs[i][j][k] >>= mesh->solution->rhs[i][j][k];
				}
			}
		}
		
		trace_off();
#else
		if(config->solver->scheme == "forward_euler"){
			equation->calc_residual(q, mesh->solution->rhs);
			update_forward_euler(q.size(), q.data(), mesh->solution->rhs.data(), mesh->solution->dt.data());
		}
		else if(config->solver->scheme == "rk4_jameson"){

			for(int order=0; order<4; order++){
				equation->calc_residual(q_tmp, mesh->solution->rhs);
				update_rk4(q.size(), q_tmp.data(), q.data(), mesh->solution->rhs.data(), mesh->solution->dt.data(), order);
			}

			set_rarray(q.size(), q.data(), q_tmp.data());
			
		}
		
		else{
			logger->critical("scheme not defined.");
		}
	
#endif
		config->profiler->update_time_residual();
		
		Tx l2normq;
		for(uint k=0; k<nq+ntrans; k++){
			l2normq = 0.0;
			for(uint i=0; i<nic; i++){
				for(uint j=0; j<njc; j++){
					l2normq += mesh->solution->rhs[i][j][k]*mesh->solution->rhs[i][j][k];
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
			break;
		}
		//		logger->debug("transport l2 {}", l2norm[4]);
		if(l2norm[0] < config->solver->tolerance && l2norm[1] < config->solver->tolerance && 0){
			logger->info("Convergence reached!");
			iomanager->write(counter);
			auto dt_main = config->profiler->current_time();
			logger->info("Final:: Step: {:08d} Time: {:.2e} Wall Time: {:.2e} CFL: {:.2e} Density Norm: {:.2e}", counter, t, dt_main, CFL, l2norm[0]);
			logger_convergence->info("{:08d} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}", counter, t, dt_main, CFL, l2norm[0], l2norm[1], l2norm[2], l2norm[3]);
			break;
		}
#if defined(ENABLE_ADOLC)
		config->profiler->reset_time_jacobian();
		Tx *q_ptr = mesh->solution->q.data();
		sparse_jac(1,mesh->solution->nt,mesh->solution->nt,mesh->solution->repeat,q_ptr,&mesh->solution->nnz,&mesh->solution->rind,&mesh->solution->cind,&mesh->solution->values,mesh->solution->options);
		config->profiler->update_time_jacobian();
		logger->debug("NNZ = {}", mesh->solution->nnz);
		
		if(counter == 0)
			linearsolver->preallocate(mesh->solution->nnz);
		Tx dt_local = 0;
		for(uint i=0; i<mesh->solution->nnz; i++){
			const auto k_idx = mesh->solution->rind[i]%(nq+ntrans);
			const auto j_idx = ((mesh->solution->rind[i]-k_idx)/(nq+ntrans))%njc;
			const auto i_idx = (mesh->solution->rind[i] - k_idx - j_idx*(nq+ntrans))/njc/(nq+ntrans);
			dt_local = mesh->solution->dt[i_idx][j_idx][k_idx];
			mesh->solution->values[i] = -mesh->solution->values[i];
			//			std::cout<<dt_local<<std::endl;
			if(mesh->solution->rind[i] == mesh->solution->cind[i]){mesh->solution->values[i] += 1.0/dt_local;}
		}
		config->profiler->reset_time_linearsolver();
		linearsolver->set_lhs(mesh->solution->nnz, mesh->solution->rind, mesh->solution->cind, mesh->solution->values);
		linearsolver->set_rhs(mesh->solution->rhs.data());
		linearsolver->solve_and_update(mesh->solution->q.data(), UNDER_RELAXATION);
		auto dt_perf = config->profiler->update_time_linearsolver();
		//q[i][j][k] = q[i][j][k] + mesh->solution->rhs[i][j][k]*dt;
		
		logger->info("Linear algebra time = {:03.2f}", dt_perf);
		
		free(mesh->solution->rind); mesh->solution->rind=nullptr;
		free(mesh->solution->cind); mesh->solution->cind=nullptr;
		free(mesh->solution->values); mesh->solution->values=nullptr;
#else
#endif
		
	
		
		t += 0;
		counter += 1;

		auto cfl_ramp =  config->solver->cfl_ramp;
		if(cfl_ramp){
			auto cfl_ramp_iteration = config->solver->cfl_ramp_iteration;
			if(counter > cfl_ramp_iteration){
				auto cfl_ramp_exponent = config->solver->cfl_ramp_exponent;
				CFL = pow(CFL, cfl_ramp_exponent);
				CFL = std::min(CFL, static_cast<Tx>(1e6));
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


		
		
		if(counter % config->io->stdout_frequency == 0){
			auto dt_main = config->profiler->current_time();
			logger->info("Step: {:08d} Time: {:.2e} Wall Time: {:.2e} CFL: {:.2e} Density Norm: {:.2e}", counter, t, dt_main, CFL, l2norm[0]);
			logger_convergence->info("{:08d} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}", counter, t, dt_main, CFL, l2norm[0], l2norm[1], l2norm[2], l2norm[3]);
		}
		
		if(counter % config->io->fileout_frequency == 0){
			iomanager->write(counter);
		}
	}
}


#endif
