#ifndef _EULEREQUATION_H
#define _EULEREQUATION_H
#include "flux.h"
#include "reconstruction.h"
#include "bc.h"
template<class T, class Tad>
class EulerEquation{
public:
	uint ni, nj, nic, njc, nq;
	Array2D<Tad> rho, u, v, p;
	Array2D<Tad> rholft_xi, ulft_xi, vlft_xi, plft_xi;
	Array2D<Tad> rhorht_xi, urht_xi, vrht_xi, prht_xi;
	Array2D<Tad> rholft_eta, ulft_eta, vlft_eta, plft_eta;
	Array2D<Tad> rhorht_eta, urht_eta, vrht_eta, prht_eta;
	Array3D<Tad> flux_xi, flux_eta;

	std::shared_ptr<Mesh<T>> mesh;
	std::shared_ptr<Config<T>> config;
	std::unique_ptr<Reconstruction<T, Tad>> reconstruction;
	std::unique_ptr<ConvectiveFlux<T, Tad>> convective_flux;
	

	std::vector<BoundaryCondition<T, Tad>*> boundaryconditions;

	void calc_residual(Array3D<Tad>& a_q, Array3D<Tad>& a_rhs);
	void calc_convective_residual(Array3D<Tad>& a_rhs);
	void calc_viscous_residual(Array3D<Tad>& a_rhs){
		auto grad_u = Array3D<Tad>(nic, njc, 3);
		auto grad_v = Array3D<Tad>(nic, njc, 3);
		auto grad_Temp = Array3D<Tad>(nic, njc, 3);
		auto Temp = Array2D<Tad>(nic+2, njc+2);
		auto Rc = 1.0;
		auto mu = 1e-4;
		auto k = mu*GAMMA/(GAMMA-1.0);
		for(int i=0; i<nic+2; i++){
			for(int j=0; j<njc+2; j++){
				Temp[i][j] = p[i][j]/rho[i][j]/Rc;
			}
		}
		mesh->calc_gradient(u, grad_u, 1, 1);
		mesh->calc_gradient(v, grad_v, 1, 1);
		mesh->calc_gradient(Temp, grad_Temp, 1, 1);

		// xi
		
		for(int i=0; i<ni; i++){
			for(int j=0; j<njc; j++){
				Tad dudx, dudy, dvdx, dvdy, ubar, vbar, dTdx, dTdy;
					auto nx = mesh->normal_chi[i][j][0];
					auto ny = mesh->normal_chi[i][j][1];
					ubar = (u[i+1][j+1] + u[i][j+1])*0.5;
					vbar = (v[i+1][j+1] + v[i][j+1])*0.5;
					if(i == 0){
					dudx = grad_u[i][j][0];
					dudy = grad_u[i][j][1];
					
					dvdx = grad_v[i][j][0];
					dvdy = grad_v[i][j][1];
					
					dTdx = grad_Temp[i][j][0];
					dTdy = grad_Temp[i][j][1];


				}

				else if(i == ni-1){
					dudx = grad_u[i-1][j][0];
					dudy = grad_u[i-1][j][1];
					
					dvdx = grad_v[i-1][j][0];
					dvdy = grad_v[i-1][j][1];
					
					dTdx = grad_Temp[i-1][j][0];
					dTdy = grad_Temp[i-1][j][1];

				}
				else{
					auto dudx_bar = (grad_u[i][j][0] + grad_u[i-1][j][0])*0.5;
					auto dudy_bar = (grad_u[i][j][1] + grad_u[i-1][j][1])*0.5;
					
					auto dvdx_bar = (grad_v[i][j][0] + grad_v[i-1][j][0])*0.5;
					auto dvdy_bar = (grad_v[i][j][1] + grad_v[i-1][j][1])*0.5;
					
					auto dTdx_bar = (grad_Temp[i][j][0] + grad_Temp[i-1][j][0])*0.5;
					auto dTdy_bar = (grad_Temp[i][j][1] + grad_Temp[i-1][j][1])*0.5;
					
					auto rijx = mesh->xc[i][j] - mesh->xc[i-1][j];
					auto rijy = mesh->yc[i][j] - mesh->yc[i-1][j];
					auto dr = std::sqrt(rijx*rijx + rijy*rijy);
					rijx = rijx/dr;
					rijy = rijy/dr;
					
					dudx = dudx_bar;// + ((u[i][j] - u[i-1][j])/dr - dudx_bar*rijx - dudy_bar*rijy)*rijx;
					dudy = dudy_bar;// + ((u[i][j] - u[i-1][j])/dr - dudx_bar*rijx - dudy_bar*rijy)*rijy;
					
					dvdx = dvdx_bar;// + ((v[i][j] - v[i-1][j])/dr - dvdx_bar*rijx - dvdy_bar*rijy)*rijx;
					dvdy = dvdy_bar;// + ((v[i][j] - v[i-1][j])/dr - dvdx_bar*rijx - dvdy_bar*rijy)*rijy;
					
					dTdx = dTdx_bar;// + ((Temp[i][j] - Temp[i-1][j])/dr - dTdx_bar*rijx - dTdy_bar*rijy)*rijx;
					dTdy = dTdy_bar;// + ((Temp[i][j] - Temp[i-1][j])/dr - dTdx_bar*rijx - dTdy_bar*rijy)*rijy;


				}

				auto tau_xy = mu*(dudy + dvdx);
				auto tau_xx = mu*(2.0*dudx - 2.0/3.0*(dudx + dvdy));
				auto tau_yy = mu*(2.0*dvdy - 2.0/3.0*(dudx + dvdy));
				auto q_x = -k*dTdx;
				auto q_y = -k*dTdy;

				flux_xi[i][j][1] = tau_xx*nx + tau_xy*ny;
				flux_xi[i][j][2] = tau_xy*nx + tau_yy*ny;
				flux_xi[i][j][3] = (ubar*tau_xx + vbar*tau_xy - q_x)*nx + (ubar*tau_xy + vbar*tau_yy - q_y)*ny;
			}
		}


				
		for(int i=0; i<nic; i++){
			for(int j=0; j<nj; j++){
				Tad dudx, dudy, dvdx, dvdy, ubar, vbar, dTdx, dTdy;
					auto nx = mesh->normal_eta[i][j][0];
					auto ny = mesh->normal_eta[i][j][1];
					ubar = (u[i+1][j+1] + u[i+1][j])*0.5;
					vbar = (v[i+1][j+1] + v[i+1][j])*0.5;

					if(j == 0){
						dudx = grad_u[i][j][0];
					dudy = grad_u[i][j][1];
					
					dvdx = grad_v[i][j][0];
					dvdy = grad_v[i][j][1];
					
					dTdx = grad_Temp[i][j][0];
					dTdy = grad_Temp[i][j][1];


				}

				else if(j == nj-1){
					dudx = grad_u[i][j-1][0];
					dudy = grad_u[i][j-1][1];
					
					dvdx = grad_v[i][j-1][0];
					dvdy = grad_v[i][j-1][1];
					
					dTdx = grad_Temp[i][j-1][0];
					dTdy = grad_Temp[i][j-1][1];


				}
				else{
					auto dudx_bar = (grad_u[i][j][0] + grad_u[i][j-1][0])*0.5;
					auto dudy_bar = (grad_u[i][j][1] + grad_u[i][j-1][1])*0.5;
					
					auto dvdx_bar = (grad_v[i][j][0] + grad_v[i][j-1][0])*0.5;
					auto dvdy_bar = (grad_v[i][j][1] + grad_v[i][j-1][1])*0.5;
					
					auto dTdx_bar = (grad_Temp[i][j][0] + grad_Temp[i][j-1][0])*0.5;
					auto dTdy_bar = (grad_Temp[i][j][1] + grad_Temp[i][j-1][1])*0.5;
					
					auto rijx = mesh->xc[i][j] - mesh->xc[i][j-1];
					auto rijy = mesh->yc[i][j] - mesh->yc[i][j-1];
					auto dr = std::sqrt(rijx*rijx + rijy*rijy);
					rijx = rijx/dr;
					rijy = rijy/dr;
					
					dudx = dudx_bar;// + ((u[i][j] - u[i][j-1])/dr - dudx_bar*rijx - dudy_bar*rijy)*rijx;
					dudy = dudy_bar;// + ((u[i][j] - u[i][j-1])/dr - dudx_bar*rijx - dudy_bar*rijy)*rijy;
					
					dvdx = dvdx_bar;// + ((v[i][j] - v[i][j-1])/dr - dvdx_bar*rijx - dvdy_bar*rijy)*rijx;
					dvdy = dvdy_bar;// + ((v[i][j] - v[i][j-1])/dr - dvdx_bar*rijx - dvdy_bar*rijy)*rijy;
					
					dTdx = dTdx_bar;// + ((Temp[i][j] - Temp[i][j-1])/dr - dTdx_bar*rijx - dTdy_bar*rijy)*rijx;
					dTdy = dTdy_bar;// + ((Temp[i][j] - Temp[i][j-1])/dr - dTdx_bar*rijx - dTdy_bar*rijy)*rijy;

				}

				auto tau_xy = mu*(dudy + dvdx);
				auto tau_xx = mu*(2.0*dudx - 2.0/3.0*(dudx + dvdy));
				auto tau_yy = mu*(2.0*dvdy - 2.0/3.0*(dudx + dvdy));
				auto q_x = -k*dTdx;
				auto q_y = -k*dTdy;

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
	void calc_source_residual(Array3D<Tad>& a_rhs){};
	void calc_primvars(Array3D<Tad>& a_q);
	void calc_boundary();
	EulerEquation(std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config<T>> val_config){
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
		
		flux_xi = Array3D<Tad>(ni, njc, 4U);
		flux_eta = Array3D<Tad>(nic, nj, 4U);

		if(config->solver->order == 1){
			reconstruction = std::make_unique<FirstOrder<T, Tad>>(ni, nj);
		}
		else if(config->solver->order == 2){
			reconstruction = std::make_unique<SecondOrder<T, Tad>>(ni, nj);
		}
		else{
			
		}

		if(config->solver->flux == "roe")
			convective_flux = std::make_unique<RoeFlux<T, Tad>>();
		else if(config->solver->flux == "ausm")
			convective_flux = std::make_unique<AUSMFlux<T, Tad>>();
		else
			spdlog::get("console")->critical("Flux not found.");
		
		for(const auto& bc : config->geometry->boundary){
			auto name = bc->name;
			auto type = bc->type;
			auto face = bc->face;
			auto start = bc->start;
			auto end = bc->end;
			uint facei;
			if(face=="left") facei = left;
			if(face=="right") facei = right;
			if(face=="bottom") facei = bottom;
			if(face=="top") facei = top;
			if(type == "freestream"){
				auto boundarycondition = new BoundaryConditionFreestream<T, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else if(type == "slipwall"){
				auto boundarycondition = new BoundaryConditionInviscidWall<T, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else if(type == "wall"){
				auto boundarycondition = new BoundaryConditionViscousWall<T, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else if(type == "wake"){
				auto boundarycondition = new BoundaryConditionWake<T, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else{
				spdlog::get("console")->info("Wrong type of BC.");
			}
		}
		
	};

	~EulerEquation(){
		for(auto&& bc : boundaryconditions)
			delete bc;
	};
};
template <class T, class Tad>
void EulerEquation<T, Tad>::calc_convective_residual(Array3D<Tad>& a_rhs){
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

template <class T, class Tad>
void EulerEquation<T, Tad>::calc_primvars(Array3D<Tad>& a_q){
	primvars<Tad>(a_q, rho, u, v, p, 1U, 1U);
}


template <class T, class Tad>
void EulerEquation<T, Tad>::calc_boundary(){
	for(auto&& bc : boundaryconditions)
		bc->apply(rho, u, v, p);
}

template <class T, class Tad>
	void EulerEquation<T, Tad>::calc_residual(Array3D<Tad>& a_q, Array3D<Tad>& a_rhs){
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
	calc_source_residual(a_rhs);
	
#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] /= mesh->volume[i][j];
			}
		}
	}

};

#endif
