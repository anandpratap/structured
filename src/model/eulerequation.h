#ifndef _EULEREQUATION_H
#define _EULEREQUATION_H
#include "common.h"
#include "flux.h"

template<class T, class Tad>
class EulerEquation{
public:
	Array2D<Tad> rho, u, v, p;
	Array2D<Tad> rholft_xi, ulft_xi, vlft_xi, plft_xi;
	Array2D<Tad> rhorht_xi, urht_xi, vrht_xi, prht_xi;
	Array2D<Tad> rholft_eta, ulft_eta, vlft_eta, plft_eta;
	Array2D<Tad> rhorht_eta, urht_eta, vrht_eta, prht_eta;
	Array3D<Tad> flux_xi, flux_eta;
	Array3D<Tad> a_q, a_rhs;
	
	std::shared_ptr<Mesh<T>> mesh;
	std::shared_ptr<Config<T>> config;
	void calc_residual(Tad *a_q_ravel, Tad *a_rhs_ravel);
 EulerEquation(std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config<T>> val_config, uint val_ni, uint val_nj):
	rho(val_ni-1+2, val_nj-1+2),
		u(val_ni-1+2, val_nj-1+2),
		v(val_ni-1+2, val_nj-1+2),
		p(val_ni-1+2, val_nj-1+2),
		rholft_xi(val_ni, val_nj-1),
		ulft_xi(val_ni, val_nj-1),
		vlft_xi(val_ni, val_nj-1),
		plft_xi(val_ni, val_nj-1),
		rhorht_xi(val_ni, val_nj-1),
		urht_xi(val_ni, val_nj-1),
		vrht_xi(val_ni, val_nj-1),
		prht_xi(val_ni, val_nj-1),

		rholft_eta(val_ni-1, val_nj),
		ulft_eta(val_ni-1, val_nj),
		vlft_eta(val_ni-1, val_nj),
		plft_eta(val_ni-1, val_nj),
		rhorht_eta(val_ni-1, val_nj),
		urht_eta(val_ni-1, val_nj),
		vrht_eta(val_ni-1, val_nj),
		prht_eta(val_ni-1, val_nj),
		flux_xi(val_ni, val_nj-1, 4U),
		flux_eta(val_ni-1, val_nj, 4U),
		a_q(val_ni-1, val_nj-1, 4U),
		a_rhs(val_ni-1, val_nj-1, 4U)
			{
		mesh = val_mesh;
		config = val_config;

		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ni = nic + 1;
		const auto nj = njc + 1;
		
	};

	~EulerEquation(){
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ni = nic + 1;
		const auto nj = njc + 1;
	};
};

template <class T, class Tad>
	void EulerEquation<T, Tad>::calc_residual(Tad *a_q_ravel, Tad *a_rhs_ravel){
	const auto ni = mesh->ni;
	const auto nj = mesh->nj;
	const auto nic = mesh->nic;
	const auto njc = mesh->njc;
	const auto nq = mesh->solution->nq;

#pragma omp parallel for
	for(uint i=0; i<nic*njc*nq; i++){
		a_q(i) = a_q_ravel[i];
	}

#pragma omp parallel for
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			//spdlog::get("console")->info("{} {}", i, j);
			Tad tmp_rho = a_q(i, j, 0);
			Tad tmp_u = a_q(i, j, 1)/tmp_rho;
			Tad tmp_v = a_q(i, j, 2)/tmp_rho;
			Tad tmp_p = (a_q(i, j, 3) - 0.5*tmp_rho*(tmp_u*tmp_u + tmp_v*tmp_v))*(GAMMA-1);
		
			rho(i+1, j+1) = tmp_rho;
			u(i+1, j+1) = tmp_u;
			v(i+1, j+1) = tmp_v;
			p(i+1, j+1) = tmp_p;
			
			
			//	primvars<Tad>(&a_q_tmp, &rho(i+1, j+1), &u(i+1, j+1), &v(i+1, j+1), &p(i+1, j+1));



			for(uint k=0; k<nq; k++){
				a_rhs(i, j, k) = 0.0;
			}
		}
	}
	//	rho.print();
	auto rho_inf =  config->freestream->rho_inf;
	auto u_inf =  config->freestream->u_inf;
	auto v_inf =  config->freestream->v_inf;
	auto p_inf =  config->freestream->p_inf;

#pragma omp parallel for
	for(uint i=0; i<nic+2; i++){
		rho(i, njc+1) = rho_inf;
		u(i, njc+1) = u_inf;
		v(i, njc+1) = v_inf;
		p(i, njc+1) = p_inf;
		
	}

#pragma omp parallel for
	for(uint j=0; j<njc+2; j++){
		rho(0, j) = rho_inf;
		u(0, j) = u_inf;
		v(0, j) = v_inf;
		p(0, j) = p_inf;
		rho(nic+1, j) = rho_inf;
		u(nic+1, j) = u_inf;
		v(nic+1, j) = v_inf;
		p(nic+1, j) = p_inf;
	}


	const auto j1 = mesh->j1;
	const auto nb = mesh->nb;

#pragma omp parallel for
	for(uint i=0; i<nb; i++){
		Tad un, ds;
		ds = mesh->normal_eta[j1-1+i][0][0]*mesh->normal_eta[j1-1+i][0][0] +
			mesh->normal_eta[j1-1+i][0][1]*mesh->normal_eta[j1-1+i][0][1];

		p(j1+i, 0) = 1.5*p(j1+i, 1) - 0.5*p(j1+i, 2);
		rho(j1+i, 0) = 1.5*rho(j1+i, 1) - 0.5*rho(j1+i, 2);
		un = u(j1+i, 1)*mesh->normal_eta[j1-1+i][0][0] + v(j1+i, 1)*mesh->normal_eta[j1-1+i][0][1];
		u(j1+i, 0) = u(j1+i, 1) - 2*un*mesh->normal_eta[j1-1+i][0][0]/ds;
		v(j1+i, 0) = v(j1+i, 1) - 2*un*mesh->normal_eta[j1-1+i][0][1]/ds;
	}

#pragma omp parallel for
	for(uint i=1; i < j1; i++){
		rho(i, 0) = rho(nic+1-i, 1);
		u(i, 0) = u(nic+1-i, 1);
		v(i, 0) = v(nic+1-i, 1);
		p(i, 0) = p(nic+1-i, 1);

		rho(nic+1-i, 0) = rho(i, 1);
		u(nic+1-i, 0) = u(i, 1);
		v(nic+1-i, 0) = v(i, 1);
		p(nic+1-i, 0) = p(i, 1);
	}


	

	if(config->solver->order == 1){
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


	if(config->solver->order == 1){
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
	Tad tmp_flux[4] = {0.0, 0.0, 0.0, 0.0};
#pragma omp parallel for
	for(uint i=0; i< ni; i++){
		for(uint j=0; j< njc; j++){
			roeflux<Tad>(mesh->normal_chi[i][j][0], mesh->normal_chi[i][j][1],
						 rholft_xi(i, j), ulft_xi(i, j), vlft_xi(i, j), plft_xi(i, j),
						 rhorht_xi(i, j), urht_xi(i, j), vrht_xi(i, j), prht_xi(i, j),
						 tmp_flux);
			
			for(uint k=0; k<4; k++){
				flux_xi(i, j, k) = tmp_flux[k];
			}
		}
	}

#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< nj; j++){
			roeflux<Tad>(mesh->normal_eta[i][j][0], mesh->normal_eta[i][j][1],
						 rholft_eta(i, j), ulft_eta(i, j), vlft_eta(i, j), plft_eta(i, j),
						 rhorht_eta(i, j), urht_eta(i, j), vrht_eta(i, j), prht_eta(i, j),
						 tmp_flux);
			for(uint k=0; k<4; k++){
				flux_eta(i, j, k) = tmp_flux[k];
			}
		}
	}

#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				//spdlog::get("console")->info("{} {} {}", i, j, k);
				//spdlog::get("console")->info("{}", a_rhs(i,j,k));
				//spdlog::get("console")->info("{}", flux_eta.get_size());
				//spdlog::get("console")->info("{}", flux_eta(i,j+1,k));
				//spdlog::get("console")->info("{}", flux_eta(i,j,k));
				a_rhs(i, j, k) += - flux_eta(i, j+1, k) + flux_eta(i, j, k);
				a_rhs(i, j, k) += -flux_xi(i+1, j, k) + flux_xi(i, j, k);
			}
		}
	}

#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs(i,j,k) /=  mesh->volume[i][j];
			}
		}
	}

#pragma omp parallel for
	for(uint i=0; i<nic*njc*nq; i++){
		a_rhs_ravel[i] = a_rhs(i);
	}
};

#endif
