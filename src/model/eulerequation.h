#ifndef _EULEREQUATION_H
#define _EULEREQUATION_H
#include "flux.h"

template<class T, class Tad>
class EulerEquation{
public:
	Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> rho, u, v, p;
	Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> rholft_xi, ulft_xi, vlft_xi, plft_xi;
	Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> rhorht_xi, urht_xi, vrht_xi, prht_xi;
	Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> rholft_eta, ulft_eta, vlft_eta, plft_eta;
	Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> rhorht_eta, urht_eta, vrht_eta, prht_eta;
	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> normal_eta_x, normal_eta_y, normal_xi_x, normal_xi_y;
	Tad ***flux_xi, ***flux_eta;
	Tad ***a_q, ***a_rhs;

	std::shared_ptr<Mesh<T>> mesh;
	std::shared_ptr<Config<T>> config;
	void calc_residual(Tad *a_q_ravel, Tad *a_rhs_ravel);
	EulerEquation(std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config<T>> val_config){
		mesh = val_mesh;
		config = val_config;

		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ni = nic + 1;
		const auto nj = njc + 1;
		
		rho = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic+2, njc+2);
		u = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic+2, njc+2);
		v = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic+2, njc+2);
		p = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic+2, njc+2);
		rholft_xi = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);
		ulft_xi = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);
		vlft_xi = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);
		plft_xi = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);

		normal_xi_x = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);
		normal_xi_y = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);
		
		rhorht_xi = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);
		urht_xi = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);
		vrht_xi = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);
		prht_xi = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(ni, njc);
		
		rholft_eta = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);
		ulft_eta = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);
		vlft_eta = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);
		plft_eta = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);
		
		rhorht_eta = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);
		urht_eta = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);
		vrht_eta = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);
		prht_eta = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);

		normal_eta_x = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);
		normal_eta_y = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>(nic, nj);
		flux_xi = allocate_3d_array<Tad>(ni, njc, 4U);
		flux_eta = allocate_3d_array<Tad>(nic, nj, 4U);
		
		a_q = allocate_3d_array<Tad>(mesh->nic, mesh->njc, mesh->solution->nq);
		a_rhs = allocate_3d_array<Tad>(mesh->nic, mesh->njc, mesh->solution->nq);

		for(int i=0; i<mesh->ni; i++){
			for(int j=0; j<mesh->njc; j++){
				normal_xi_x(i,j) = mesh->normal_chi[i][j][0];
				normal_xi_y(i,j) = mesh->normal_chi[i][j][1];
			}
		}

		for(int i=0; i<mesh->nic; i++){
			for(int j=0; j<mesh->nj; j++){
				normal_eta_x(i,j) = mesh->normal_eta[i][j][0];
				normal_eta_y(i,j) = mesh->normal_eta[i][j][1];
			}
		}
		
	};

	~EulerEquation(){
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ni = nic + 1;
		const auto nj = njc + 1;
		release_3d_array(flux_xi, ni, njc, 4U);
		release_3d_array(flux_eta, nic, nj, 4U);
		release_3d_array(a_q, mesh->nic, mesh->njc, mesh->solution->nq);
		release_3d_array(a_rhs, mesh->nic, mesh->njc, mesh->solution->nq);
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
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			for(uint k=0; k<nq; k++){
				a_q[i][j][k] = a_q_ravel[i*njc*nq + j*nq + k];
			}
		}
	}
	
	Tad ***q = a_q;

#pragma omp parallel for
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			//spdlog::get("console")->info("{} {}", i, j);
			primvars<Tad>(q[i][j], &rho(i+1, j+1), &u(i+1, j+1), &v(i+1, j+1), &p(i+1, j+1));
			for(uint k=0; k<nq; k++){
				a_rhs[i][j][k] = 0.0;
			}
		}
	}
	auto rho_inf =  config->freestream->rho_inf;
	auto u_inf =  config->freestream->u_inf;
	auto v_inf =  config->freestream->v_inf;
	auto p_inf =  config->freestream->p_inf;

	rho.block(0, njc+1, nic+2, 1).setConstant(rho_inf);
	u.block(0, njc+1, nic+2, 1).setConstant(u_inf);
	v.block(0, njc+1, nic+2, 1).setConstant(v_inf);
	p.block(0, njc+1, nic+2, 1).setConstant(p_inf);

	rho.block(0, 0, 1, njc+2).setConstant(rho_inf);
	u.block(0, 0, 1, njc+2).setConstant(u_inf);
	v.block(0, 0, 1, njc+2).setConstant(v_inf);
	p.block(0, 0, 1, njc+2).setConstant(p_inf);

	rho.block(nic+1, 0, 1, njc+2).setConstant(rho_inf);
	u.block(nic+1, 0, 1, njc+2).setConstant(u_inf);
	v.block(nic+1, 0, 1, njc+2).setConstant(v_inf);
	p.block(nic+1, 0, 1, njc+2).setConstant(p_inf);

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

#pragma omp parallel for
	for(uint i=0; i< ni; i++){
		for(uint j=0; j< njc; j++){
			roeflux(normal_xi_x(i, j), normal_xi_y(i, j),
					rholft_xi(i, j), ulft_xi(i, j), vlft_xi(i, j), plft_xi(i, j),
					rhorht_xi(i, j), urht_xi(i, j), vrht_xi(i, j), prht_xi(i, j),
					flux_xi[i][j]);
		}
	}
	
#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< nj; j++){
			roeflux(normal_eta_x(i, j), normal_eta_y(i, j),
					rholft_eta(i, j), ulft_eta(i, j), vlft_eta(i, j), plft_eta(i, j),
					rhorht_eta(i, j), urht_eta(i, j), vrht_eta(i, j), prht_eta(i, j),
					flux_eta[i][j]);
		}
	}

#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] -= (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
				a_rhs[i][j][k] -= (flux_xi[i+1][j][k] - flux_xi[i][j][k]);
			}
		}
	}

#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] /= mesh->volume[i][j];
			}
		}
	}

#pragma omp parallel for
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			for(uint k=0; k<nq; k++){
				a_rhs_ravel[i*njc*nq + j*nq + k] = a_rhs[i][j][k];
			}
		}
	}
};

#endif
