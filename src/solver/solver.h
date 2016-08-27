#ifndef _SOLVER_H
#define _SOLVER_H
#include "common.h"
#include "utils.h"
#include "mesh.h"
#include "linearsolver.h"
#include "config.h"
#include "io.h"
#include "eulerequation.h"

template<class T>
void set_rarray(int size, T* __restrict__ dest, T* __restrict__ src){
#pragma omp parallel for
	for(auto i=0; i<size; i++){
		dest[i] = src[i];
	}
}

template<class T>
void update_forward_euler(int size, T* __restrict__ q, T* __restrict__ rhs, T* __restrict__ dt){
#pragma omp parallel for
	for(auto i=0; i<size; i++){
		q[i] = q[i] + rhs[i]*dt[i];
	}
}

template<class T, class To>
	void update_rk4(int size, T* __restrict__ q_i, T* __restrict__ q, T* __restrict__ rhs, T* __restrict__ dt, To order){
#pragma omp parallel for
	for(auto i=0; i<size; i++){
		q_i[i] = q[i] + rhs[i]*dt[i]/(4.0 - order);
	}
}

template<class T, class Tad>
class Solver{
public:
	int nnz;
	int repeat = 0;
	unsigned int *rind = nullptr;
	unsigned int *cind = nullptr;
	double *values = nullptr;
	int options[4] = {0,0,0,0};
	Timer timer_la;
	Timer timer_main;
	Timer timer_residual;
	std::shared_ptr<Mesh<T>> mesh;
	void solve();
	T UNDER_RELAXATION;
	T CFL;
	std::shared_ptr<Config<T>> config;
	std::string label;
	std::shared_ptr<spdlog::logger> logger;
	std::shared_ptr<spdlog::logger> logger_convergence;	
	std::shared_ptr<IOManager<T>> iomanager;

	std::shared_ptr<EulerEquation<T, Tad>> equation;

#if defined(ENABLE_ARMA)
	std::shared_ptr<LinearSolverArma<T>> linearsolver;
#endif
#if defined(ENABLE_EIGEN)
	std::shared_ptr<LinearSolverEigen<T>> linearsolver;
#endif
#if defined(ENABLE_PETSC)
	std::shared_ptr<LinearSolverPetsc<T>> linearsolver;
#endif

	uint nt;
	Array3D<T> rhs;
	T **lhs;
	Array3D<T> dt;
	Array3D<T> q;
	Array3D<T> q_tmp;
	Array3D<Tad> a_q, a_rhs;
	
	Solver(std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config<T>> config);
	~Solver();
	void copy_from_solution();
	void copy_to_solution();
	void calc_dt();
	void initialize();
};
template <class T, class Tad>
void Solver<T, Tad>::calc_dt(){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;

#pragma omp parallel for
	for(uint i=0; i<mesh->nic; i++){
		for(uint j=0; j<mesh->njc; j++){
			T rho = q[i][j][0];
			T u = q[i][j][1]/rho;
			T v = q[i][j][2]/rho;
			T rhoE = q[i][j][3];
			T p = (rhoE - 0.5*rho*(u*u + v*v))*(GAMMA-1.0);
			T lambda = std::sqrt(GAMMA*p/rho) + std::abs(u) + std::abs(v);
			for(int k=0; k<nq; k++)
				dt[i][j][k] = std::min(mesh->ds_eta[i][j], mesh->ds_chi[i][j])/lambda*CFL;
		}
	}

}
template <class T, class Tad>
void Solver<T, Tad>::copy_from_solution(){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;

	set_rarray(q.size(), q.data(), mesh->solution->q.data());
}

template <class T, class Tad>
void Solver<T, Tad>::copy_to_solution(){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	set_rarray(q.size(), mesh->solution->q.data(), q.data());
}

template <class T, class Tad>
void Solver<T, Tad>::initialize(){
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
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			for(uint k=0; k<nq; k++){
				q_tmp[i][j][k] = mesh->solution->q[i][j][k];
			}
		}
	}
	
}

template <class T, class Tad>
Solver<T, Tad>::Solver(std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config<T>> val_config){
	config = val_config;
	timer_la = Timer();
	timer_main = Timer();
	timer_residual = Timer();
	mesh = val_mesh;
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	nt = nic*njc*nq;
	
	dt = Array3D<T>(nic, njc, nq);
	rhs = Array3D<T>(nic, njc, nq);
	q = Array3D<T>(nic, njc, nq);
	q_tmp = Array3D<T>(nic, njc, nq);
	a_q = Array3D<Tad>(nic, njc, nq);
	a_rhs = Array3D<Tad>(nic, njc, nq);

#if defined(ENABLE_ARMA)
	linearsolver = std::make_shared<LinearSolverArma<T>>(mesh, config);
#endif

#if defined(ENABLE_EIGEN)
	linearsolver = std::make_shared<LinearSolverEigen<T>>(mesh, config);
#endif

#if defined(ENABLE_PETSC)
	linearsolver = std::make_shared<LinearSolverPetsc<T>>(mesh, config);
#endif


	logger_convergence = spdlog::basic_logger_mt("convergence", "history.dat", true);
	logger_convergence->info(" ");
	logger = spdlog::get("console");

	CFL = config->solver->cfl;
	UNDER_RELAXATION = config->solver->under_relaxation;
	label = config->io->label;


	equation = std::make_shared<EulerEquation<T, Tad>>(mesh, config);
	iomanager = std::make_shared<IOManager<T>>(mesh, config);
}
template <class T, class Tad>
Solver<T, Tad>::~Solver(){
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	
}


template <class T, class Tad>
void Solver<T, Tad>::solve(){
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	T l2norm[nq] = {1e10};

	uint counter = 0;
	T t = 0.0;
	initialize();
	copy_from_solution();
	logger->info("Welcome to structured!");


	config->profiler->timer_main->reset();

	while(1){
		calc_dt();
		timer_residual.reset();
		config->profiler->reset_time_residual();

#if defined(ENABLE_ADOLC)

		
		trace_on(1);
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq; k++){
					a_q[i][j][k] <<= q[i][j][k];
				}
			}
		}
		equation->calc_residual(a_q, a_rhs);
		
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq; k++){
					a_rhs[i][j][k] >>= rhs[i][j][k];
				}
			}
		}
		
		trace_off();
#else
		if(config->solver->scheme == "forward_euler"){
			equation->calc_residual(q, rhs);
			update_forward_euler(q.size(), q.data(), rhs.data(), dt.data());
		}
		else if(config->solver->scheme == "rk4_jameson"){

			for(int order=0; order<4; order++){
				equation->calc_residual(q_tmp, rhs);
				update_rk4(q.size(), q_tmp.data(), q.data(), rhs.data(), dt.data(), order);
			}

			set_rarray(q.size(), q.data(), q_tmp.data());
			
		}
		
		else{
			logger->critical("scheme not defined.");
		}
	
#endif
		config->profiler->update_time_residual();
		float dt_perfs = timer_residual.diff();

		T l2normq;
		for(uint k=0; k<nq; k++){
			l2normq = 0.0;
			for(uint i=0; i<nic; i++){
				for(uint j=0; j<njc; j++){
					l2normq += rhs[i][j][k]*rhs[i][j][k];
				}
			}
			l2norm[k] = sqrt(l2normq);
		}
		
		if(counter > config->solver->iteration_max){
			logger->info("Max iteration reached!");
			break;
		}
		if(l2norm[0] < config->solver->tolerance){
			logger->info("Convergence reached!");
			copy_to_solution();
			iomanager->write(counter);
			float dt_main = timer_main.diff();
			logger->info("Final:: Step: {:08d} Time: {:.2e} Wall Time: {:.2e} CFL: {:.2e} Density Norm: {:.2e}", counter, t, dt_main, CFL, l2norm[0]);
			break;
		}
#if defined(ENABLE_ADOLC)
		config->profiler->reset_time_jacobian();
		T *q_ptr = q.data();
		sparse_jac(1,nt,nt,repeat,q_ptr,&nnz,&rind,&cind,&values,options);
		config->profiler->update_time_jacobian();
		logger->debug("NNZ = {}", nnz);
		
		if(counter == 0)
			linearsolver->preallocate(nnz);
		T dt_local = 0;
		for(uint i=0; i<nnz; i++){
			const auto k_idx = rind[i]%nq;
			const auto j_idx = ((rind[i]-k_idx)/nq)%njc;
			const auto i_idx = (rind[i] - k_idx - j_idx*nq)/njc/nq;
			dt_local = dt[i_idx][j_idx][k_idx];
			values[i] = -values[i];
			//			std::cout<<dt_local<<std::endl;
			if(rind[i] == cind[i]){values[i] += 1.0/dt_local;}
		}
		timer_la.reset();
		config->profiler->reset_time_linearsolver();
		linearsolver->set_jac(nnz, rind, cind, values);
		linearsolver->set_rhs(rhs.data());
		linearsolver->solve_and_update(q.data(), UNDER_RELAXATION);
		config->profiler->update_time_linearsolver();

		//q[i][j][k] = q[i][j][k] + rhs[i][j][k]*dt;
		
		float dt_perf = timer_la.diff();
		logger->info("Linear algebra time = {:03.2f}", dt_perf);
		
		free(rind); rind=nullptr;
		free(cind); cind=nullptr;
		free(values); values=nullptr;
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
				CFL = std::min(CFL, static_cast<T>(1e6));
			}
		}

		auto under_relaxation_ramp =  config->solver->under_relaxation_ramp;
		if(under_relaxation_ramp){
			auto under_relaxation_ramp_iteration = config->solver->under_relaxation_ramp_iteration;
			if(counter > under_relaxation_ramp_iteration){
				auto under_relaxation_ramp_exponent = config->solver->under_relaxation_ramp_exponent;
				UNDER_RELAXATION = pow(UNDER_RELAXATION, under_relaxation_ramp_exponent);
				UNDER_RELAXATION = std::min(UNDER_RELAXATION, static_cast<T>(10.0));
			}
		}


		
		
		if(counter % config->io->stdout_frequency == 0){
			float dt_main = timer_main.diff();
			logger->info("Step: {:08d} Time: {:.2e} Wall Time: {:.2e} CFL: {:.2e} Density Norm: {:.2e}", counter, t, dt_main, CFL, l2norm[0]);
			logger_convergence->info("{:08d} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}", counter, t, dt_main, CFL, l2norm[0], l2norm[1], l2norm[2], l2norm[3]);
		}
		
		if(counter % config->io->fileout_frequency == 0){
			copy_to_solution();
			iomanager->write(counter);
		}
	}
}


#endif
