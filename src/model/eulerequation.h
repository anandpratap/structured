#ifndef _EULEREQUATION_H
#define _EULEREQUATION_H
#include "flux.h"
#include "reconstruction.h"
#include "bc.h"
#include "common.h"

template<class Tx, class Tad>
class EulerEquation{
public:
	uint ni, nj, nic, njc, nq;
	Array2D<Tad> rho, u, v, p, T;
	Array3D<Tad> grad_u, grad_v, grad_T;
	Array2D<Tad> rholft_xi, ulft_xi, vlft_xi, plft_xi;
	Array2D<Tad> rhorht_xi, urht_xi, vrht_xi, prht_xi;
	Array2D<Tad> rholft_eta, ulft_eta, vlft_eta, plft_eta;
	Array2D<Tad> rhorht_eta, urht_eta, vrht_eta, prht_eta;
	Array3D<Tad> flux_xi, flux_eta;

	Array3D<Tad> grad_u_xi, grad_u_eta;
	Array3D<Tad> grad_v_xi, grad_v_eta;
	Array3D<Tad> grad_T_xi, grad_T_eta;
	Array2D<Tad> u_bar_xi, u_bar_eta;
	Array2D<Tad> v_bar_xi, v_bar_eta;
	
	std::shared_ptr<Mesh<Tx>> mesh;
	std::shared_ptr<Config<Tx>> config;
	std::unique_ptr<Reconstruction<Tx, Tad>> reconstruction;
	std::unique_ptr<ConvectiveFlux<Tx, Tad>> convective_flux;
	
	std::unique_ptr<BoundaryContainer<Tx, Tad>> boundary_container;
   
	//std::unique_ptr<TransportEquation<Tx, Tad>> transport;	
	void calc_residual(Array3D<Tad>& a_q, Array3D<Tad>& a_rhs);
	void calc_convective_residual(Array3D<Tad>& a_rhs);
	void calc_viscous_residual(Array3D<Tad>& a_rhs){
		auto mu = config->freestream->mu_inf;
		auto p_inf = config->freestream->p_inf;
		auto pr_inf = config->freestream->pr_inf;
		auto rho_inf = config->freestream->rho_inf;
		auto T_inf = config->freestream->T_inf;
		auto Rc = p_inf/rho_inf/T_inf;
		auto k = mu*GAMMA*Rc/(GAMMA-1.0)/pr_inf;
	
		mesh->calc_gradient(u, grad_u_xi, grad_u_eta);
		mesh->calc_gradient(v, grad_v_xi, grad_v_eta);
		mesh->calc_gradient(T, grad_T_xi, grad_T_eta);
		mesh->calc_face(u, u_bar_xi, u_bar_eta);
		mesh->calc_face(v, v_bar_xi, v_bar_eta);
		// xi
		
		for(int i=0; i<ni; i++){
			for(int j=0; j<njc; j++){
				const Tx nx = mesh->normal_chi[i][j][0];
				const Tx ny = mesh->normal_chi[i][j][1];
				const Tad dudx = grad_u_xi[i][j][0];
				const Tad dudy = grad_u_xi[i][j][1];

				const Tad dvdx = grad_v_xi[i][j][0];
				const Tad dvdy = grad_v_xi[i][j][1];

				const Tad dTdx = grad_T_xi[i][j][0];
				const Tad dTdy = grad_T_xi[i][j][1];

				const Tad ubar = u_bar_xi[i][j];
				const Tad vbar = v_bar_xi[i][j];

				const Tad tau_xy = mu*(dudy + dvdx);
				const Tad tau_xx = mu*(2.0*dudx - 2.0/3.0*(dudx + dvdy));
				const Tad tau_yy = mu*(2.0*dvdy - 2.0/3.0*(dudx + dvdy));
				const Tad q_x = -k*dTdx;
				const Tad q_y = -k*dTdy;

				flux_xi[i][j][0] = 0.0;
				flux_xi[i][j][1] = tau_xx*nx + tau_xy*ny;
				flux_xi[i][j][2] = tau_xy*nx + tau_yy*ny;
				flux_xi[i][j][3] = (ubar*tau_xx + vbar*tau_xy - q_x)*nx + (ubar*tau_xy + vbar*tau_yy - q_y)*ny;
			}
		}

		for(int i=0; i<nic; i++){
			for(int j=0; j<nj; j++){
				const Tx nx = mesh->normal_eta[i][j][0];
				const Tx ny = mesh->normal_eta[i][j][1];
				const Tad dudx = grad_u_eta[i][j][0];
				const Tad dudy = grad_u_eta[i][j][1];

				const Tad dvdx = grad_v_eta[i][j][0];
				const Tad dvdy = grad_v_eta[i][j][1];

				const Tad dTdx = grad_T_eta[i][j][0];
				const Tad dTdy = grad_T_eta[i][j][1];

				const Tad ubar = u_bar_eta[i][j];
				const Tad vbar = v_bar_eta[i][j];

				const Tad tau_xy = mu*(dudy + dvdx);
				const Tad tau_xx = mu*(2.0*dudx - 2.0/3.0*(dudx + dvdy));
				const Tad tau_yy = mu*(2.0*dvdy - 2.0/3.0*(dudx + dvdy));
				const Tad q_x = -k*dTdx;
				const Tad q_y = -k*dTdy;

				flux_eta[i][j][0] = 0.0;
				flux_eta[i][j][1] = tau_xx*nx + tau_xy*ny;
				flux_eta[i][j][2] = tau_xy*nx + tau_yy*ny;
				flux_eta[i][j][3] = (ubar*tau_xx + vbar*tau_xy - q_x)*nx + (ubar*tau_xy + vbar*tau_yy - q_y)*ny;
			}
		}

		
		for(uint i=0; i< nic; i++){
			for(uint j=0; j< njc; j++){
				for(uint k=1; k<mesh->solution->nq; k++){
					a_rhs[i][j][k] += (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
					a_rhs[i][j][k] += (flux_xi[i+1][j][k] - flux_xi[i][j][k]);
				}
			}
		}
	};
	void calc_source_residual(Array3D<Tad>& a_q, Array3D<Tad>& a_rhs){
		for(uint i=0; i< nic; i++){
			for(uint j=0; j< njc; j++){
				a_rhs[i][j][1] += -config->solver->dpdx*mesh->volume[i][j];
				a_rhs[i][j][2] += -config->solver->dpdy*mesh->volume[i][j];
			}
		}
	};
	void calc_primvars(Array3D<Tad>& a_q);
	void calc_boundary();
	EulerEquation(std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
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
		
		rholft_xi = Array2D<Tad>(ni, njc);
		ulft_xi = Array2D<Tad>(ni, njc);
		vlft_xi = Array2D<Tad>(ni, njc);
		plft_xi = Array2D<Tad>(ni, njc);
		
		rhorht_xi = Array2D<Tad>(ni, njc);
		urht_xi = Array2D<Tad>(ni, njc);
		vrht_xi = Array2D<Tad>(ni, njc);
		prht_xi = Array2D<Tad>(ni, njc);
		
		rholft_eta = Array2D<Tad>(nic, nj);
		ulft_eta = Array2D<Tad>(nic, nj);
		vlft_eta = Array2D<Tad>(nic, nj);
		plft_eta = Array2D<Tad>(nic, nj);
		
		rhorht_eta = Array2D<Tad>(nic, nj);
		urht_eta = Array2D<Tad>(nic, nj);
		vrht_eta = Array2D<Tad>(nic, nj);
		prht_eta = Array2D<Tad>(nic, nj);

		grad_u_xi = Array3D<Tad>(ni, njc, 2);
		grad_u_eta = Array3D<Tad>(nic, nj, 2);

		grad_v_xi = Array3D<Tad>(ni, njc, 2);
		grad_v_eta = Array3D<Tad>(nic, nj, 2);

		grad_T_xi = Array3D<Tad>(ni, njc, 2);
		grad_T_eta = Array3D<Tad>(nic, nj, 2);

		u_bar_xi = Array2D<Tad>(ni, njc);
		v_bar_xi = Array2D<Tad>(ni, njc);

		u_bar_eta = Array2D<Tad>(nic, nj);
		v_bar_eta = Array2D<Tad>(nic, nj);
			
		grad_u = Array3D<Tad>(nic, njc, 3);
		grad_v = Array3D<Tad>(nic, njc, 3);
		grad_T = Array3D<Tad>(nic, njc, 3);
		

		flux_xi = Array3D<Tad>(ni, njc, 4U);
		flux_eta = Array3D<Tad>(nic, nj, 4U);

		//transport = std::make_unique<TransportEquation<Tx, Tad>>(mesh, config);

		
		if(config->solver->order == 1){
			reconstruction = std::make_unique<FirstOrder<Tx, Tad>>(ni, nj);
		}
		else if(config->solver->order == 2){
			reconstruction = std::make_unique<SecondOrder<Tx, Tad>>(ni, nj);
		}
		else{
			
		}

		boundary_container = std::make_unique<BoundaryContainer<Tx, Tad>>(config->filename, mesh, config);
		
		if(config->solver->flux == "roe")
			convective_flux = std::make_unique<RoeFlux<Tx, Tad>>();
		else if(config->solver->flux == "ausm")
			convective_flux = std::make_unique<AUSMFlux<Tx, Tad>>();
		else
			spdlog::get("console")->critical("Flux not found.");
				
	};

	~EulerEquation(){
	};
};
template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_convective_residual(Array3D<Tad>& a_rhs){
	convective_flux->evaluate(mesh->normal_chi,
							  rholft_xi, ulft_xi, vlft_xi, plft_xi,
							  rhorht_xi, urht_xi, vrht_xi, prht_xi,
							  flux_xi);
	convective_flux->evaluate(mesh->normal_eta,
							  rholft_eta, ulft_eta, vlft_eta, plft_eta,
							  rhorht_eta, urht_eta, vrht_eta, prht_eta,
							  flux_eta);
#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] -= (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
				a_rhs[i][j][k] -= (flux_xi[i+1][j][k] - flux_xi[i][j][k]);
			}
		}
	}
}

template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_primvars(Array3D<Tad>& a_q){
	primvars<Tad>(a_q, rho, u, v, p, T, 1U, 1U);
}


template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_boundary(){
	boundary_container->apply(rho, u, v, p, T);
}

template <class Tx, class Tad>
	void EulerEquation<Tx, Tad>::calc_residual(Array3D<Tad>& a_q, Array3D<Tad>& a_rhs){
	a_rhs.fill(0.0);

	calc_primvars(a_q);
	calc_boundary();
	
	reconstruction->evaluate_xi(rho, rholft_xi, rhorht_xi);
	reconstruction->evaluate_xi(u, ulft_xi, urht_xi);
	reconstruction->evaluate_xi(v, vlft_xi, vrht_xi);
	reconstruction->evaluate_xi(p, plft_xi, prht_xi);

	reconstruction->evaluate_eta(rho, rholft_eta, rhorht_eta);
	reconstruction->evaluate_eta(u, ulft_eta, urht_eta);
	reconstruction->evaluate_eta(v, vlft_eta, vrht_eta);
	reconstruction->evaluate_eta(p, plft_eta, prht_eta);
	
	calc_convective_residual(a_rhs);
	calc_viscous_residual(a_rhs);
	calc_source_residual(a_q, a_rhs);
	//transport->calc_residual(a_q, a_rhs, u, v,
	//						 ulft_xi, urht_xi,
	//						 ulft_eta, urht_eta,
	//						 vlft_xi, vrht_xi,
	//						 vlft_eta, vrht_eta);
#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq+mesh->solution->ntrans; k++){
				a_rhs[i][j][k] /= mesh->volume[i][j];
			}
		}
	}

};

#endif
