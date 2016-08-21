#ifndef _SOLVER_H
#define _SOLVER_H
#include "common.h"
#include "utils.h"
#include "mesh.h"
#include "linearsolver.h"

template<class T>
class Solver{
public:
	int nnz;
	int repeat = 0;
	unsigned int *rind = nullptr;
	unsigned int *cind = nullptr;
	double *values = nullptr;
	int options[4] = {0,0,0,1};
	Timer timer;
	Timer timer_main;
	Mesh<T> *mesh;
	void calc_residual();
	void solve();
	adouble ***a_q, ***a_rhs;
	double UNDER_RELAXATION;
	double CFL;
	std::string label;
	std::shared_ptr<cpptoml::table> config;
	std::shared_ptr<spdlog::logger> logger;
	std::shared_ptr<spdlog::logger> logger_convergence;	
#if defined(ENABLE_ARMA)
	LinearSolverArma *linearsolver;
#endif
#if defined(ENABLE_EIGEN)
	LinearSolverEigen *linearsolver;
#endif
#if defined(ENABLE_PETSC)
	LinearSolverPetsc *linearsolver;
#endif
	
	T *rhs;
	T **lhs;
	T *dt;
	T *q;
	adouble *a_q_ravel, *a_rhs_ravel;

	Solver(Mesh<T> *val_mesh);
	~Solver();
	void copy_from_solution();
	void copy_to_solution();
	void calc_dt();
	void initialize();
};
template <class T>
void Solver<T>::calc_dt(){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;

	for(uint i=0; i<mesh->nic; i++){
		for(uint j=0; j<mesh->njc; j++){
			double rho = q[i*njc*nq + j*nq + 0];
			double u = q[i*njc*nq + j*nq + 1]/rho;
			double v = q[i*njc*nq + j*nq + 2]/rho;
			double rhoE = q[i*njc*nq + j*nq + 3];
			double p = (rhoE - 0.5*rho*(u*u + v*v))*(GAMMA-1.0);
			double lambda = sqrt(GAMMA*p/rho) + abs(u) + abs(v);
			dt[i*njc + j] = std::min(mesh->ds_eta[i][j], mesh->ds_chi[i][j])/lambda*CFL;
		}
	}

}
template <class T>
void Solver<T>::copy_from_solution(){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;

	for(uint i=0; i<mesh->nic; i++){
		for(uint j=0; j<mesh->njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				q[i*njc*nq + j*nq + k] = mesh->solution->q[i][j][k];
			}
		}
	}
}

template <class T>
void Solver<T>::copy_to_solution(){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;

	for(uint i=0; i<mesh->nic; i++){
		for(uint j=0; j<mesh->njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				mesh->solution->q[i][j][k] = q[i*njc*nq + j*nq + k];
			}
		}
	}
}

template <class T>
void Solver<T>::initialize(){
	uint nic = mesh->nic;
	uint njc = mesh->njc;

	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			auto rho_inf =  config->get_qualified_as<double>("freestream.rho_inf").value_or(1.0);
			auto u_inf =  config->get_qualified_as<double>("freestream.u_inf").value_or(0.0);
			auto v_inf =  config->get_qualified_as<double>("freestream.v_inf").value_or(0.0);
			auto p_inf =  config->get_qualified_as<double>("freestream.p_inf").value_or(1.0/1.4);
			mesh->solution->q[i][j][0] = rho_inf;
			mesh->solution->q[i][j][1] = rho_inf*u_inf;
			mesh->solution->q[i][j][2] = rho_inf*v_inf;
			mesh->solution->q[i][j][3] = p_inf/(GAMMA-1.0) + 0.5*rho_inf*(u_inf*u_inf + v_inf*v_inf);
		}
	}
	if(config->get_qualified_as<bool>("io.restart").value_or(false)){
		std::string filename_out = label + ".out";
		read_restart_file(mesh, filename_out);
	}
}

template <class T>
Solver<T>::Solver(Mesh<T> *val_mesh){
	timer = Timer();
	timer_main = Timer();
	mesh = val_mesh;
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
 
	a_q = allocate_3d_array<adouble>(mesh->nic, mesh->njc, mesh->solution->nq);
	a_rhs = allocate_3d_array<adouble>(mesh->nic, mesh->njc, mesh->solution->nq);
	
	dt = allocate_1d_array<T>(nic*njc);
	rhs = allocate_1d_array<T>(nic*njc*nq);
	q = allocate_1d_array<T>(nic*njc*nq);
	a_q_ravel = allocate_1d_array<adouble>(nic*njc*nq);
	a_rhs_ravel = allocate_1d_array<adouble>(nic*njc*nq);

#if defined(ENABLE_ARMA)
	linearsolver = new LinearSolverArma(mesh);
#endif

#if defined(ENABLE_EIGEN)
	linearsolver = new LinearSolverEigen(mesh);
#endif

#if defined(ENABLE_PETSC)
	linearsolver = new LinearSolverPetsc(mesh);
#endif


	spdlog::set_pattern("%v");
	logger_convergence = spdlog::basic_logger_mt("convergence", "history.dat", true);
	logger_convergence->info(" ");
	logger = spdlog::stdout_logger_mt("console", true);
	logger->set_level(spdlog::level::debug);
	config = cpptoml::parse_file("config.inp");
	CFL = config->get_qualified_as<double>("solver.cfl").value_or(1.0);
	UNDER_RELAXATION = config->get_qualified_as<double>("solver.under_relaxation").value_or(1.0);
	label = config->get_qualified_as<std::string>("io.label").value_or("flow");
}
template <class T>
Solver<T>::~Solver(){
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;

	release_3d_array(a_q, mesh->nic, mesh->njc, mesh->solution->nq);
	release_3d_array(a_rhs, mesh->nic, mesh->njc, mesh->solution->nq);
	release_1d_array(q, nic*njc*nq);
	release_1d_array(a_q_ravel, nic*njc*nq);
	release_1d_array(a_rhs_ravel, nic*njc*nq);
	release_1d_array(rhs, nic*njc*nq);
	release_1d_array(dt, nic*njc);
	delete linearsolver;
}


template <class T>
void Solver<T>::solve(){
#if defined(ENABLE_EIGEN)
	Eigen::SparseLU<Eigen::SparseMatrix<double>> eigen_solver;
	//Eigen::SuperLU<Eigen::SparseMatrix<double>> eigen_solver;
#endif
	
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	double l2norm[nq] = {1e10};

	uint counter = 0;
	double t = 0.0;
	initialize();
	copy_from_solution();
	logger->info("Welcome to structured!");
	while(1){
		trace_on(1);
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq; k++){
					a_q_ravel[i*njc*nq + j*nq + k] <<= q[i*njc*nq + j*nq + k];
				}
			}
		}
		
		calc_residual();

		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq; k++){
					a_rhs_ravel[i*njc*nq + j*nq + k] >>= rhs[i*njc*nq + j*nq + k];
				}
			}
		}
		
		trace_off();


		for(int j=0; j<nq; j++){
			l2norm[j] = 0.0;
			for(int i=0; i<nic*njc; i++){
				l2norm[j] += rhs[i*nq+j]*rhs[i*nq+j];
			}
			l2norm[j] = sqrt(l2norm[j]);
		}
		if(l2norm[0] < 1e-8){
			logger->info("Convergence reached!");
			float dt_main = timer_main.diff();
			logger->info("Final:: Step: {:08d} Time: {:.2e} Wall Time: {:.2e} CFL: {:.2e} Density Norm: {:.2e}", counter, t, dt_main, CFL, l2norm[0]);
			break;
		}


		sparse_jac(1,nic*njc*nq,nic*njc*nq,repeat,q,&nnz,&rind,&cind,&values,options);
		logger->debug("NNZ = {}", nnz);

		if(counter == 0)
			linearsolver->preallocate(nnz);
		calc_dt();
		double dt_local = 0;
		for(uint i=0; i<nnz; i++){
			dt_local = dt[rind[i]/nq];
			values[i] = -values[i];
			//			std::cout<<dt_local<<std::endl;
			if(rind[i] == cind[i]){values[i] += 1.0/dt_local;}
		}
		timer.reset();
		linearsolver->set_jac(nnz, rind, cind, values);
		linearsolver->set_rhs(rhs);
		linearsolver->solve_and_update(q, UNDER_RELAXATION);
		
	//q[i][j][k] = q[i][j][k] + rhs[i][j][k]*dt;
		
		float dt_perf = timer.diff();
		logger->info("Linear algebra time = {:03.2f}", dt_perf);

		free(rind); rind=nullptr;
		free(cind); cind=nullptr;
		free(values); values=nullptr;
		

		/* for(uint i=0; i<nic; i++){ */
		/* 	for(uint j=0; j<njc; j++){ */
		/* 		for(uint k=0; k<nq; k++){ */
		/* 			//q[i][j][k] = q[i][j][k] + rhs[i][j][k]*dt; */
		/* 			q[i][j][k] = q[i][j][k] + arma_dq(i*njc*nq + j*nq + k, 0); */
		/* 		} */
		/* 	} */
		/* } */
		t += 0;
		counter += 1;

		auto cfl_ramp =  config->get_qualified_as<bool>("solver.cfl_ramp").value_or(false);
		if(cfl_ramp){
			auto cfl_ramp_iteration = config->get_qualified_as<int64_t>("solver.cfl_ramp_iteration").value_or(20);
			if(counter > cfl_ramp_iteration){
				auto cfl_ramp_exponent = config->get_qualified_as<double>("solver.cfl_ramp_exponent").value_or(1.1);
				CFL = pow(CFL, cfl_ramp_exponent);
				CFL = std::min(CFL, 1e6);
			}
		}

		auto under_relaxation_ramp =  config->get_qualified_as<bool>("solver.under_relaxation_ramp").value_or(false);
		if(under_relaxation_ramp){
			auto under_relaxation_ramp_iteration = config->get_qualified_as<int64_t>("solver.under_relaxation_ramp_iteration").value_or(20);
			if(counter > under_relaxation_ramp_iteration){
				auto under_relaxation_ramp_exponent = config->get_qualified_as<double>("solver.under_relaxation_ramp_exponent").value_or(1.1);
				UNDER_RELAXATION = pow(UNDER_RELAXATION, under_relaxation_ramp_exponent);
				UNDER_RELAXATION = std::min(UNDER_RELAXATION, 10.0);
			}
		}

		
		if(counter % 1 == 0){
			float dt_main = timer_main.diff();
			logger->info("Step: {:08d} Time: {:.2e} Wall Time: {:.2e} CFL: {:.2e} Density Norm: {:.2e}", counter, t, dt_main, CFL, l2norm[0]);
			logger_convergence->info("{:08d} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}", counter, t, dt_main, CFL, l2norm[0], l2norm[1], l2norm[2], l2norm[3]);
			copy_to_solution();

			std::string filename_tec = label + ".tec";
			std::string filename_npy = label + ".npy";
			std::string filename_out = label + ".out";
			std::string filename_surface = label + ".surface";
			write_solution(mesh, filename_tec);
			write_solution_npy(mesh, filename_npy);
			write_restart_file(mesh, filename_out);
			write_surface_file(mesh, filename_surface);
		}		
	}
}

template <class T>
void Solver<T>::calc_residual(){
	//std::cout<<"calc_res"<<std::endl;
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			for(uint k=0; k<nq; k++){
				a_q[i][j][k] = a_q_ravel[i*njc*nq + j*nq + k];
			}
		}
	}
	
	adouble ***q = a_q;
	static adouble **rho = allocate_2d_array<adouble>(nic+2, njc+2);
	static adouble **u = allocate_2d_array<adouble>(nic+2, njc+2);
	static adouble **v = allocate_2d_array<adouble>(nic+2, njc+2);
	static adouble **p = allocate_2d_array<adouble>(nic+2, njc+2);
	
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			primvars<adouble>(q[i][j], &rho[i+1][j+1], &u[i+1][j+1], &v[i+1][j+1], &p[i+1][j+1]);
			// std::cout<<rho[i][j]<<std::endl;
			// std::cout<<u[i][j]<<std::endl;
			// std::cout<<v[i][j]<<std::endl;
			// std::cout<<p[i][j]<<std::endl;
			for(uint k=0; k<nq; k++){
				a_rhs[i][j][k] = 0.0;
			}
		}
	}
	auto rho_inf =  config->get_qualified_as<double>("freestream.rho_inf").value_or(1.0);
	auto u_inf =  config->get_qualified_as<double>("freestream.u_inf").value_or(0.0);
	auto v_inf =  config->get_qualified_as<double>("freestream.v_inf").value_or(0.0);
	auto p_inf =  config->get_qualified_as<double>("freestream.p_inf").value_or(1.0/1.4);
	
	for(uint i=0; i<nic+2; i++){
		rho[i][njc+1] = rho_inf;
		u[i][njc+1] = u_inf;
		v[i][njc+1] = v_inf;
		p[i][njc+1] = p_inf;
		
	}

	for(uint j=0; j<njc+2; j++){
		rho[0][j] = rho_inf;
		u[0][j] = u_inf;
		v[0][j] = v_inf;
		p[0][j] = p_inf;
		rho[nic+1][j] = rho_inf;
		u[nic+1][j] = u_inf;
		v[nic+1][j] = v_inf;
		p[nic+1][j] = p_inf;
	}


	uint j1 = mesh->j1;
	uint nb = mesh->nb;

	adouble un, ds;

	for(uint i=0; i<nb; i++){
		ds = mesh->normal_eta[j1-1+i][0][0]*mesh->normal_eta[j1-1+i][0][0] +
			mesh->normal_eta[j1-1+i][0][1]*mesh->normal_eta[j1-1+i][0][1];

		p[j1+i][0] = 1.5*p[j1+i][1] - 0.5*p[j1+i][2];
		rho[j1+i][0] = 1.5*rho[j1+i][1] - 0.5*rho[j1+i][2];
		un = u[j1+i][1]*mesh->normal_eta[j1-1+i][0][0] + v[j1+i][1]*mesh->normal_eta[j1-1+i][0][1];
		u[j1+i][0] = u[j1+i][1] - 2*un*mesh->normal_eta[j1-1+i][0][0]/ds;
		v[j1+i][0] = v[j1+i][1] - 2*un*mesh->normal_eta[j1-1+i][0][1]/ds;
	}

	for(uint i=1; i < j1; i++){
		rho[i][0] = rho[nic+1-i][1];
		u[i][0] = u[nic+1-i][1];
		v[i][0] = v[nic+1-i][1];
		p[i][0] = p[nic+1-i][1];

		rho[nic+1-i][0] = rho[i][1];
		u[nic+1-i][0] = u[i][1];
		v[nic+1-i][0] = v[i][1];
		p[nic+1-i][0] = p[i][1];
	}


	
	static adouble** rholft_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** ulft_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** vlft_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** plft_xi = allocate_2d_array<adouble>(ni, njc);

	static adouble** rhorht_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** urht_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** vrht_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** prht_xi = allocate_2d_array<adouble>(ni, njc);

	if(config->get_qualified_as<int64_t>("solver.order").value_or(1) == 1){
		first_order_xi(ni, nj, rho, rholft_xi, rhorht_xi);
		first_order_xi(ni, nj, u, ulft_xi, urht_xi);
		first_order_xi(ni, nj, v, vlft_xi, vrht_xi);
		first_order_xi(ni, nj, p, plft_xi, prht_xi);
	}
	else{
		second_order_xi(ni, nj, rho, rholft_xi, rhorht_xi);
		second_order_xi(ni, nj, u, ulft_xi, urht_xi);
		second_order_xi(ni, nj, v, vlft_xi, vrht_xi);
		second_order_xi(ni, nj, p, plft_xi, prht_xi);
	
	}

	static adouble** rholft_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** ulft_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** vlft_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** plft_eta = allocate_2d_array<adouble>(nic, nj);

	static adouble** rhorht_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** urht_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** vrht_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** prht_eta = allocate_2d_array<adouble>(nic, nj);

	if(config->get_qualified_as<int64_t>("solver.order").value_or(1) == 1){
	first_order_eta(ni, nj, rho, rholft_eta, rhorht_eta);
	first_order_eta(ni, nj, u, ulft_eta, urht_eta);
	first_order_eta(ni, nj, v, vlft_eta, vrht_eta);
	first_order_eta(ni, nj, p, plft_eta, prht_eta);
	}else{
	second_order_eta(ni, nj, rho, rholft_eta, rhorht_eta);
	second_order_eta(ni, nj, u, ulft_eta, urht_eta);
	second_order_eta(ni, nj, v, vlft_eta, vrht_eta);
	second_order_eta(ni, nj, p, plft_eta, prht_eta);

	}

	static adouble*** flux_xi = allocate_3d_array<adouble>(ni, njc, 4U);
	static adouble*** flux_eta = allocate_3d_array<adouble>(nic, nj, 4U);
	
	for(uint i=0; i< ni; i++){
		for(uint j=0; j< njc; j++){
			roeflux<adouble>(mesh->normal_chi[i][j][0], mesh->normal_chi[i][j][1],
					   rholft_xi[i][j], ulft_xi[i][j], vlft_xi[i][j], plft_xi[i][j],
					   rhorht_xi[i][j], urht_xi[i][j], vrht_xi[i][j], prht_xi[i][j],
					   flux_xi[i][j]);
			//std::cout<<flux_xi[i][j][0]<<" "<<rholft_xi[i][j]<<" "<<rhorht_xi[i][j]<<std::endl;
			//std::cout<<flux_xi[i][j][1]<<" "<<ulft_xi[i][j]<<" "<<urht_xi[i][j]<<std::endl;
			//std::cout<<flux_xi[i][j][2]<<" "<<vlft_xi[i][j]<<" "<<vrht_xi[i][j]<<std::endl;
			//std::cout<<flux_xi[i][j][3]<<" "<<plft_xi[i][j]<<" "<<prht_xi[i][j]<<std::endl;

		}
	}

	for(uint i=0; i< nic; i++){
		for(uint j=0; j< nj; j++){
			roeflux<adouble>(mesh->normal_eta[i][j][0], mesh->normal_eta[i][j][1],
					   rholft_eta[i][j], ulft_eta[i][j], vlft_eta[i][j], plft_eta[i][j],
					   rhorht_eta[i][j], urht_eta[i][j], vrht_eta[i][j], prht_eta[i][j],
					   flux_eta[i][j]);
		}
	}

	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] -= (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
				a_rhs[i][j][k] -= (flux_xi[i+1][j][k] - flux_xi[i][j][k]);
			}
		}
	}

	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] /= mesh->volume[i][j];
			}
		}
	}

	
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			for(uint k=0; k<nq; k++){
				a_rhs_ravel[i*njc*nq + j*nq + k] = a_rhs[i][j][k];
			}
		}
	}
	
	
};

#endif
