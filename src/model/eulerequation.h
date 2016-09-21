#ifndef _EULEREQUATION_H
#define _EULEREQUATION_H
#include "flux.h"
#include "reconstruction.h"
#include "bc.h"
#include "common.h"
#include "fluid.h"
#include "def_eulerequation.h"
template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_viscous_residual(Array3D<Tad>& a_rhs){
	diffusive_flux->evaluate(mesh->normal_chi.const_ref(), grad_u_chi.const_ref(), grad_v_chi.const_ref(), grad_T_chi.const_ref(),
							 u_bar_chi.const_ref(), v_bar_chi.const_ref(), mu_bar_chi.const_ref(), k_bar_chi.const_ref(), flux_chi);
	diffusive_flux->evaluate(mesh->normal_eta.const_ref(), grad_u_eta.const_ref(), grad_v_eta.const_ref(), grad_T_eta.const_ref(),
							 u_bar_eta.const_ref(), v_bar_eta.const_ref(), mu_bar_eta.const_ref(), k_bar_eta.const_ref(), flux_eta);
	
	for(size_t i=0; i< nic; i++){
		for(size_t j=0; j< njc; j++){
			for(size_t k=1; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] += (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
				a_rhs[i][j][k] += (flux_chi[i+1][j][k] - flux_chi[i][j][k]);
			}
		}
	}
	
};
template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_source_residual(const Array3D<const Tad>& a_q, Array3D<Tad>& a_rhs){
	for(size_t i=0; i< nic; i++){
		for(size_t j=0; j< njc; j++){
			a_rhs[i][j][1] += -config->solver->dpdx*mesh->volume[i][j];
			a_rhs[i][j][2] += -config->solver->dpdy*mesh->volume[i][j];
		}
	}
};
template <class Tx, class Tad>
EulerEquation<Tx, Tad>::EulerEquation(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
	mesh = val_mesh;
	config = val_config;

	ni = mesh->ni;
	nj = mesh->nj;
	nq = mesh->solution->nq;
	nic = ni - 1;
	njc = nj - 1;
		
	rho = Array2D<Tad>(nic+2, njc+2);
	u = Array2D<Tad>(nic+2, njc+2);
	v = Array2D<Tad>(nic+2, njc+2);
	p = Array2D<Tad>(nic+2, njc+2);
	T = Array2D<Tad>(nic+2, njc+2);
	mu = Array2D<Tad>(nic+2, njc+2);
	k = Array2D<Tad>(nic+2, njc+2);
		
	rholft_chi = Array2D<Tad>(ni, njc);
	ulft_chi = Array2D<Tad>(ni, njc);
	vlft_chi = Array2D<Tad>(ni, njc);
	plft_chi = Array2D<Tad>(ni, njc);
		
	rhorht_chi = Array2D<Tad>(ni, njc);
	urht_chi = Array2D<Tad>(ni, njc);
	vrht_chi = Array2D<Tad>(ni, njc);
	prht_chi = Array2D<Tad>(ni, njc);
		
	rholft_eta = Array2D<Tad>(nic, nj);
	ulft_eta = Array2D<Tad>(nic, nj);
	vlft_eta = Array2D<Tad>(nic, nj);
	plft_eta = Array2D<Tad>(nic, nj);
		
	rhorht_eta = Array2D<Tad>(nic, nj);
	urht_eta = Array2D<Tad>(nic, nj);
	vrht_eta = Array2D<Tad>(nic, nj);
	prht_eta = Array2D<Tad>(nic, nj);

	grad_u_chi = Array3D<Tad>(ni, njc, 2);
	grad_u_eta = Array3D<Tad>(nic, nj, 2);

	grad_v_chi = Array3D<Tad>(ni, njc, 2);
	grad_v_eta = Array3D<Tad>(nic, nj, 2);

	grad_T_chi = Array3D<Tad>(ni, njc, 2);
	grad_T_eta = Array3D<Tad>(nic, nj, 2);

	mu_bar_chi = Array2D<Tad>(ni, njc);
	k_bar_chi = Array2D<Tad>(ni, njc);
	u_bar_chi = Array2D<Tad>(ni, njc);
	v_bar_chi = Array2D<Tad>(ni, njc);
	T_bar_chi = Array2D<Tad>(ni, njc);

	mu_bar_eta = Array2D<Tad>(nic, nj);
	k_bar_eta = Array2D<Tad>(nic, nj);
	u_bar_eta = Array2D<Tad>(nic, nj);
	v_bar_eta = Array2D<Tad>(nic, nj);
	T_bar_eta = Array2D<Tad>(nic, nj);
			
	grad_u = Array3D<Tad>(nic, njc, 3);
	grad_v = Array3D<Tad>(nic, njc, 3);
	grad_T = Array3D<Tad>(nic, njc, 3);
		

	flux_chi = Array3D<Tad>(ni, njc, 4U);
	flux_eta = Array3D<Tad>(nic, nj, 4U);

	//transport = std::make_unique<TransportEquation<Tx, Tad>>(mesh, config);

		
	if(config->solver->order == 1){
		reconstruction_rhs = std::make_unique<ReconstructionFirstOrder<Tx, Tad>>(ni, nj);
	}
	else if(config->solver->order == 2){
		reconstruction_rhs = std::make_unique<ReconstructionSecondOrder<Tx, Tad>>(ni, nj);
	}
	else{
		spdlog::get("console")->critical("Reconstruction not found.");
	}

	if(config->solver->lhs_order == 1){
		reconstruction_lhs = std::make_unique<ReconstructionFirstOrder<Tx, Tad>>(ni, nj);
	}
	else if(config->solver->lhs_order == 2){
		reconstruction_lhs = std::make_unique<ReconstructionSecondOrder<Tx, Tad>>(ni, nj);
	}
	else{
		spdlog::get("console")->critical("Reconstruction not found.");
	}

		
		
	if(config->solver->flux == "roe")
		convective_flux = std::make_unique<ConvectiveFluxRoe<Tx, Tad>>();
	else if(config->solver->flux == "ausm")
		convective_flux = std::make_unique<ConvectiveFluxAUSM<Tx, Tad>>();
	else
		spdlog::get("console")->critical("Flux not found.");


	diffusive_flux = std::make_unique<DiffusiveFluxGreenGauss<Tx, Tad>>();
	boundary_container = std::make_unique<BoundaryContainer<Tx, Tad>>(config->filename, mesh, config, mesh->fluid_model);
};

template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_convective_residual(Array3D<Tad>& a_rhs){
	convective_flux->evaluate(mesh->normal_chi.const_ref(),
							  rholft_chi.const_ref(), ulft_chi.const_ref(), vlft_chi.const_ref(), plft_chi.const_ref(),
							  rhorht_chi.const_ref(), urht_chi.const_ref(), vrht_chi.const_ref(), prht_chi.const_ref(),
							  flux_chi);
	convective_flux->evaluate(mesh->normal_eta.const_ref(),
							  rholft_eta.const_ref(), ulft_eta.const_ref(), vlft_eta.const_ref(), plft_eta.const_ref(),
							  rhorht_eta.const_ref(), urht_eta.const_ref(), vrht_eta.const_ref(), prht_eta.const_ref(),
							  flux_eta);
#pragma omp parallel for
	for(size_t i=0; i< nic; i++){
		for(size_t j=0; j< njc; j++){
			for(size_t k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] -= (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
				a_rhs[i][j][k] -= (flux_chi[i+1][j][k] - flux_chi[i][j][k]);
			}
		}
	}
}

template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_primvars(const Array3D<const Tad>& a_q){
	mesh->fluid_model->primvars(a_q.const_ref(), rho, u, v, p, T, 1U, 1U);
}


template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_boundary(){
	boundary_container->apply(rho, u, v, p, T);
}

template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_intermediates(const Array3D<const Tad>& a_q){
	calc_primvars(a_q.const_ref());
	calc_boundary();
	
	reconstruction->evaluate_chi(rho.const_ref(), rholft_chi, rhorht_chi);
	reconstruction->evaluate_chi(u.const_ref(), ulft_chi, urht_chi);
	reconstruction->evaluate_chi(v.const_ref(), vlft_chi, vrht_chi);
	reconstruction->evaluate_chi(p.const_ref(), plft_chi, prht_chi);

	reconstruction->evaluate_eta(rho.const_ref(), rholft_eta, rhorht_eta);
	reconstruction->evaluate_eta(u.const_ref(), ulft_eta, urht_eta);
	reconstruction->evaluate_eta(v.const_ref(), vlft_eta, vrht_eta);
	reconstruction->evaluate_eta(p.const_ref(), plft_eta, prht_eta);

	if(config->freestream->if_viscous){
		for(size_t i=0; i<nic+2; i++){
			for(size_t j=0; j<njc+2; j++){
				mu[i][j] = mesh->fluid_model->get_laminar_viscosity(T[i][j]);
				k[i][j] = mesh->fluid_model->get_thermal_conductivity(T[i][j]);
			}
		}
		
		mesh->calc_gradient(u.const_ref(), grad_u_chi, grad_u_eta);
		mesh->calc_gradient(v.const_ref(), grad_v_chi, grad_v_eta);
		mesh->calc_gradient(T.const_ref(), grad_T_chi, grad_T_eta);
		mesh->calc_face(u.const_ref(), u_bar_chi, u_bar_eta);
		mesh->calc_face(v.const_ref(), v_bar_chi, v_bar_eta);
		mesh->calc_face(T.const_ref(), T_bar_chi, T_bar_eta);
		mesh->calc_face(mu.const_ref(), mu_bar_chi, mu_bar_eta);
		mesh->calc_face(k.const_ref(), k_bar_chi, k_bar_eta);
	}
};

template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_residual(const Array3D<const Tad>& a_q, Array3D<Tad>& a_rhs, const bool lhs){
	if(lhs){
		reconstruction = reconstruction_lhs;
	}
	else{
		reconstruction = reconstruction_rhs;
	}
	
	a_rhs.fill(0.0);

	calc_intermediates(a_q.const_ref());

	calc_convective_residual(a_rhs);

	if(config->freestream->if_viscous){
		calc_viscous_residual(a_rhs);
	}

	calc_source_residual(a_q.const_ref(), a_rhs);

	// divide by volume
#pragma omp parallel for
	for(size_t i=0; i< nic; i++){
		for(size_t j=0; j< njc; j++){
			for(size_t k=0; k<mesh->solution->nq+mesh->solution->ntrans; k++){
				a_rhs[i][j][k] /= mesh->volume[i][j];
			}
		}
	}

};



template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_dt(const Tx cfl){
	size_t nic = mesh->nic;
	size_t njc = mesh->njc;
	size_t nq = mesh->solution->nq;

#pragma omp parallel for
	for(size_t i=0; i<mesh->nic; i++){
		for(size_t j=0; j<mesh->njc; j++){
			Tx rho = mesh->solution->q[i][j][0];
			Tx u = mesh->solution->q[i][j][1]/rho;
			Tx v = mesh->solution->q[i][j][2]/rho;
			Tx rhoE = mesh->solution->q[i][j][3];
			Tx p = (rhoE - 0.5*rho*(u*u + v*v))*(GAMMA-1.0);
			Tx lambda = sqrt(GAMMA*p/rho) + fabs(u) + fabs(v);
			Tx len_min = std::min(mesh->ds_eta[i][j], mesh->ds_chi[i][j]);
			Tx mu = config->freestream->mu_inf;
			for(size_t k=0; k<nq; k++)
				mesh->solution->dt[i][j][k] = cfl/(lambda/len_min + 2.0*mu/len_min/len_min);
			//dt[i][j][4] = dt[i][j][0];
		}
	}

}

template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::initialize(){
	size_t nic = mesh->nic;
	size_t njc = mesh->njc;

	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
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
		mesh->iomanager->read_restart();
	}
	size_t nq = mesh->solution->nq;
	size_t ntrans = mesh->solution->ntrans;
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			for(size_t k=0; k<nq+ntrans; k++){
				mesh->solution->q_tmp[i][j][k] = mesh->solution->q[i][j][k];
			}
		}
	}
	
}

#endif
