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
	void calc_viscous_residual(Array3D<Tad>& a_rhs){};
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

		convective_flux = std::make_unique<RoeFlux<T, Tad>>();
		
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
