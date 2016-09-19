#ifndef _MESH_H
#define _MESH_H
#include "common.h"
#include "solution.h"
#include "config.h"
#include "eulerequation.h"
#include "io.h"
#include <memory>
#include "linearsolver.h"

template<class Tx>
class Solution;

template <class Tx>
class Mesh: public std::enable_shared_from_this<Mesh<Tx>>{
 public:
	std::string label;
	std::shared_ptr<Config<Tx>> config;
	uint ni, nj;
	uint nic, njc;
	uint j1, nb;
	Array2D<Tx> xv, yv;
	Array2D<Tx> xc, yc;

	Array3D<Tx> normal_eta;
	Array3D<Tx> normal_chi;

	Array2D<Tx> ds_eta;
	Array2D<Tx> ds_chi;

	Array2D<Tx> volume;
	Array2D<Tx> chi_x, chi_y, eta_x, eta_y;
	Array2D<Tx> x_chi, y_chi, x_eta, y_eta;
	
	std::shared_ptr<Solution<Tx>> solution;
	std::shared_ptr<IOManager<Tx>> iomanager;
	std::shared_ptr<EulerEquation<Tx, adouble>> equation;
	std::shared_ptr<FluidModel<Tx, adouble>> fluid_model;

#if defined(ENABLE_ARMA)
	std::shared_ptr<LinearSolverArma<Tx>> linearsolver;
#endif
#if defined(ENABLE_EIGEN)
	std::shared_ptr<LinearSolverEigen<Tx>> linearsolver;
#endif
#if defined(ENABLE_PETSC)
	std::shared_ptr<LinearSolverPetsc<Tx>> linearsolver;
#endif

public:
	Mesh(std::shared_ptr<Config<Tx>> config);

	Mesh(std::shared_ptr<Mesh<Tx>> mesh, const uint nskipi=0, const uint nskipj=0, const uint refine=0);
	~Mesh();

	void calc_metrics();

	template<class Tad>
	void calc_gradient(Array2D<Tad>& q, Array3D<Tad> &grad_q, uint skipi=0, uint skipj=0);

	template<class Tad>
	void calc_gradient(Array2D<Tad>& q, Array3D<Tad>& grad_chi, Array3D<Tad>& grad_eta);

	template<class Tad>
	void calc_face(Array2D<Tad>& q, Array2D<Tad>& q_chi, Array2D<Tad>& q_eta);	
	void setup();

	void simple_loader(std::string filename);
	void plot3d_loader(std::string filename);
};



template<class Tx>
template<class Tad>
void Mesh<Tx>::calc_face(Array2D<Tad>& q, Array2D<Tad>& q_chi, Array2D<Tad>& q_eta){
	for(int i=0; i<ni; i++){
		for(int j=0; j<njc; j++){
			const uint b = 1;
			const Tad qleft = q[i+b-1][j+b];
			const Tad qright = q[i+b][j+b];
			const Tad qtop = 0.25*(qleft+qright+q[i+b-1][j+b+1]+q[i+b][j+b+1]);
			const Tad qbottom = 0.25*(qleft+qright+q[i+b-1][j+b-1]+q[i+b][j+b-1]);
			q_chi[i][j] = 0.25*(qleft+qright+qtop+qbottom);
		}
	}

	for(int i=0; i<nic; i++){
		for(int j=0; j<nj; j++){
			const uint b = 1;
			const Tad qtop = q[i+b][j+b];
			const Tad qbottom = q[i+b][j+b-1];
			const Tad qleft = 0.25*(qtop + qbottom + q[i+b-1][j+b] + q[i+b-1][j+b-1]);
			const Tad qright = 0.25*(qtop + qbottom + q[i+b+1][j+b] + q[i+b+1][j+b-1]);
			q_eta[i][j] = 0.25*(qleft+qright+qtop+qbottom);
		}
	}
}

template<class Tx>
template<class Tad>
void Mesh<Tx>::calc_gradient(Array2D<Tad>& q, Array3D<Tad>& grad_chi, Array3D<Tad>& grad_eta){
	// chi
	Tx normal_top[2];
	Tx normal_bottom[2];
	Tx normal_left[2];
	Tx normal_right[2];
	Tx volume_;
	
	for(int i=0; i<ni; i++){
		for(int j=0; j<njc; j++){
			auto nx = normal_chi[i][j][0];
			auto ny = normal_chi[i][j][1];
			const uint b = 1;
			const Tad qleft = q[i+b-1][j+b];
			const Tad qright = q[i+b][j+b];
			const Tad qtop = 0.25*(qleft+qright+q[i+b-1][j+b+1]+q[i+b][j+b+1]);
			const Tad qbottom = 0.25*(qleft+qright+q[i+b-1][j+b-1]+q[i+b][j+b-1]);
			const Tad qbar = 0.25*(qleft+qright+qtop+qbottom);
			if(i == 0){

				volume_ = volume[i][j];
				for(int k=0; k<2; k++){
					normal_top[k] = normal_eta[i][j+1][k];
					normal_bottom[k] = normal_eta[i][j][k];
					normal_right[k] = 0.5*(normal_chi[i][j][k] + normal_chi[i+1][j][k]);
					normal_left[k] = normal_chi[i][j][k];
				}
			}
			
			else if(i == ni-1){
				volume_ = volume[i-1][j];
				for(int k=0; k<2; k++){
					normal_top[k] = normal_eta[i-1][j+1][k];
					normal_bottom[k] = normal_eta[i-1][j][k];
					normal_right[k] = normal_chi[i][j][k];
					normal_left[k] = 0.5*(normal_chi[i][j][k] + normal_chi[i-1][j][k]);
				}
			}
			else{
				volume_ = 0.5*(volume[i][j] + volume[i-1][j]);
				for(int k=0; k<2; k++){
					normal_top[k] = 0.5*(normal_eta[i][j+1][k] + normal_eta[i-1][j+1][k]);
					normal_bottom[k] = 0.5*(normal_eta[i][j][k] + normal_eta[i-1][j][k]);
					normal_right[k] = 0.5*(normal_chi[i][j][k] + normal_chi[i+1][j][k]);
					normal_left[k] = 0.5*(normal_chi[i][j][k] + normal_chi[i-1][j][k]);
				}
			}
			grad_chi[i][j][0] = (normal_top[0]*qtop - normal_bottom[0]*qbottom + normal_right[0]*qright - normal_left[0]*qleft)/volume_;
			grad_chi[i][j][1] = (normal_top[1]*qtop - normal_bottom[1]*qbottom + normal_right[1]*qright - normal_left[1]*qleft)/volume_;
		}
	}

	for(int i=0; i<nic; i++){
			for(int j=0; j<nj; j++){
				auto nx = normal_eta[i][j][0];
				auto ny = normal_eta[i][j][1];
				const uint b = 1;
				const Tad qtop = q[i+b][j+b];
				const Tad qbottom = q[i+b][j+b-1];
				
				const Tad qleft = 0.25*(qtop + qbottom + q[i+b-1][j+b] + q[i+b-1][j+b-1]);
				const Tad qright = 0.25*(qtop + qbottom + q[i+b+1][j+b] + q[i+b+1][j+b-1]);
				const Tad qbar = 0.25*(qleft+qright+qtop+qbottom);
				if(j == 0){
					volume_ = volume[i][j];
					for(int k=0; k<2; k++){
						normal_top[k] = 0.5*(normal_eta[i][j][k] + normal_eta[i][j+1][k]);
						normal_bottom[k] = normal_eta[i][j][k];
						normal_left[k] = normal_chi[i][j][k];
						normal_right[k] = normal_chi[i+1][j][k];
					}
				}
				else if(j == nj-1){
					for(int k=0; k<2; k++){
						normal_top[k] = normal_eta[i][j][k];
						normal_bottom[k] = 0.5*(normal_eta[i][j][k] + normal_eta[i][j-1][k]);
						normal_left[k] = normal_chi[i][j-1][k];
						normal_right[k] = normal_chi[i+1][j-1][k];
					}
					volume_ = volume[i][j-1];
				}
				else{
					for(int k=0; k<2; k++){
							normal_top[k] = 0.5*(normal_eta[i][j][k] + normal_eta[i][j+1][k]);
							normal_bottom[k] = 0.5*(normal_eta[i][j][k] + normal_eta[i][j-1][k]);
							normal_left[k] = 0.5*(normal_chi[i][j][k] + normal_chi[i][j-1][k]);
							normal_right[k] = 0.5*(normal_chi[i+1][j][k] + normal_chi[i+1][j-1][k]);
					}
					volume_ = 0.5*(volume[i][j] + volume[i][j-1]);
					
				}
			grad_eta[i][j][0] = (normal_top[0]*qtop - normal_bottom[0]*qbottom + normal_right[0]*qright - normal_left[0]*qleft)/volume_;
			grad_eta[i][j][1] = (normal_top[1]*qtop - normal_bottom[1]*qbottom + normal_right[1]*qright - normal_left[1]*qleft)/volume_;	
			}
	}
}

template<class Tx>
void Mesh<Tx>::simple_loader(std::string filename){
	std::ifstream infile(filename);
	infile >> std::fixed >> std::setprecision(20);
	for(uint j=0; j<nj; j++){
		for(uint i=0; i<ni; i++){  
			infile >> xv[i][j] >> yv[i][j];
		}
	}
	infile.close();
}

template<class Tx>
void Mesh<Tx>::plot3d_loader(std::string filename){
	std::ifstream infile(filename);
	infile >> std::fixed >> std::setprecision(20);
	int nblock, ni_, nj_;
	infile >> nblock;
	infile >> ni_ >> nj_;
	assert(nblock == 1);
	spdlog::get("console")->debug("{} {} {} {}", ni, ni_, nj, nj_);
	assert(ni == ni_);
	assert(nj == nj_);

	for(uint j=0; j<nj; j++){
		for(uint i=0; i<ni; i++){  
			infile >> xv[i][j];
		}
	}
	
	for(uint j=0; j<nj; j++){
		for(uint i=0; i<ni; i++){  
			infile >> yv[i][j];
		}
	}
	infile.close();
}

template<class Tx>
void Mesh<Tx>::calc_metrics(){
	auto reta = Array3D<Tx>(nic, nj, 2U);
	auto rchi = Array3D<Tx>(ni, njc, 2U);

	for(uint i=0; i<nic; i++){
		for(uint j=0; j<nj; j++){
			reta[i][j][0] = xv[i+1][j] - xv[i][j];
			reta[i][j][1] = yv[i+1][j] - yv[i][j];
			normal_eta[i][j][0] = -reta[i][j][1];
			normal_eta[i][j][1] = reta[i][j][0];
		}
	}

	for(uint i=0; i<ni; i++){
		for(uint j=0; j<njc; j++){
			rchi[i][j][0] = xv[i][j+1] - xv[i][j];
			rchi[i][j][1] = yv[i][j+1] - yv[i][j];
			normal_chi[i][j][0] = rchi[i][j][1];
			normal_chi[i][j][1] = -rchi[i][j][0];
		}
	}

	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			volume[i][j] = 0.5*(reta[i][j][0]*rchi[i][j][1] - rchi[i][j][0]*reta[i][j][1]
								+ reta[i][j+1][0]*rchi[i+1][j][1] - rchi[i+1][j][0]*reta[i][j+1][1]);
			
			xc[i][j] = 0.25*(xv[i][j] + xv[i+1][j] + xv[i][j+1] + xv[i+1][j+1]);
			yc[i][j] = 0.25*(yv[i][j] + yv[i+1][j] + yv[i][j+1] + yv[i+1][j+1]);

			ds_eta[i][j] = sqrt(std::pow(normal_eta[i][j][0],2) + std::pow(normal_eta[i][j][1],2));
			ds_chi[i][j] = sqrt(std::pow(normal_chi[i][j][0],2) + std::pow(normal_chi[i][j][1],2));
		}
	}

	for(uint i=1; i<nic-1; i++){
		for(uint j=1; j<njc-1; j++){
			x_chi[i][j] = (xc[i][j+1] - xc[i][j-1])/2.0; 
			y_chi[i][j] = (yc[i][j+1] - yc[i][j-1])/2.0;

			x_eta[i][j] = (xc[i+1][j] - xc[i-1][j])/2.0;
			y_eta[i][j] = (yc[i+1][j] - yc[i-1][j])/2.0;
		}
	}
	
	for(uint i=1; i<nic-1; i++){
		uint j = 0;
		x_chi[i][j] = (xc[i][j+1] - xc[i][j]);
		y_chi[i][j] = (yc[i][j+1] - yc[i][j]);
		x_eta[i][j] = (xc[i+1][j] - xc[i-1][j])/2.0;
		y_eta[i][j] = (yc[i+1][j] - yc[i-1][j])/2.0;
		
		j = njc-1;
		x_chi[i][j] = (xc[i][j] - xc[i][j-1]);
		y_chi[i][j] = (yc[i][j] - yc[i][j-1]);
		x_eta[i][j] = (xc[i+1][j] - xc[i-1][j])/2.0;
		y_eta[i][j] = (yc[i+1][j] - yc[i-1][j])/2.0;
	}

	for(uint j=1; j<njc-1; j++){
		uint i = 0;
		x_chi[i][j] = (xc[i][j+1] - xc[i][j-1])/2.0; 
		y_chi[i][j] = (yc[i][j+1] - yc[i][j-1])/2.0;
		x_eta[i][j] = (xc[i+1][j] - xc[i][j]);
		y_eta[i][j] = (yc[i+1][j] - yc[i][j]);

		i = nic-1;
		x_chi[i][j] = (xc[i][j+1] - xc[i][j-1])/2.0; 
		y_chi[i][j] = (yc[i][j+1] - yc[i][j-1])/2.0;
		x_eta[i][j] = (xc[i][j] - xc[i-1][j]);
		y_eta[i][j] = (yc[i][j] - yc[i-1][j]);
	}

	{
		uint i=0; uint j=0;
		x_chi[i][j] = (xc[i][j+1] - xc[i][j]);
		y_chi[i][j] = (yc[i][j+1] - yc[i][j]);
		x_eta[i][j] = (xc[i+1][j] - xc[i][j]);
		y_eta[i][j] = (yc[i+1][j] - yc[i][j]);

		i=nic-1; j=0;
		x_chi[i][j] = (xc[i][j+1] - xc[i][j]);
		y_chi[i][j] = (yc[i][j+1] - yc[i][j]);
		x_eta[i][j] = (xc[i][j] - xc[i-1][j]);
		y_eta[i][j] = (yc[i][j] - yc[i-1][j]);

		i=nic-1; j=njc-1;
		x_chi[i][j] = (xc[i][j] - xc[i][j-1]);
		y_chi[i][j] = (yc[i][j] - yc[i][j-1]);
		x_eta[i][j] = (xc[i][j] - xc[i-1][j]);
		y_eta[i][j] = (yc[i][j] - yc[i-1][j]);

		i=0; j=njc-1;
		x_chi[i][j] = (xc[i][j] - xc[i][j-1]);
		y_chi[i][j] = (yc[i][j] - yc[i][j-1]);
		x_eta[i][j] = (xc[i+1][j] - xc[i][j]);
		y_eta[i][j] = (yc[i+1][j] - yc[i][j]);
	}
	for(int i=0; i<nic; i++){
		for(int j=0; j<njc; j++){
			Tx ijac = 1.0/(x_chi[i][j]*y_eta[i][j] - x_eta[i][j]*y_chi[i][j]); 
			chi_x[i][j] = y_eta[i][j]*ijac;
			eta_x[i][j] = -y_chi[i][j]*ijac;
			chi_y[i][j] = -x_eta[i][j]*ijac;
			eta_y[i][j] = x_chi[i][j]*ijac;
		}
	}
};

template<class Tx>
template<class Tad>
void Mesh<Tx>::calc_gradient(Array2D<Tad>& q, Array3D<Tad> &grad_q, uint skipi, uint skipj){
	Tad q_chi, q_eta;
	for(uint i=1; i<nic-1; i++){
		for(uint j=1; j<njc-1; j++){
			q_chi = (q[i+skipi][j+skipj+1] - q[i+skipi][j+skipj-1])/2.0; 
			q_eta = (q[i+skipi+1][j+skipj] - q[i+skipi-1][j+skipj])/2.0;
			grad_q[i][j][0] = chi_x[i][j]*q_chi + eta_x[i][j]*q_eta;
			grad_q[i][j][1] = chi_y[i][j]*q_chi + eta_y[i][j]*q_eta;
		}
	}
	
	for(uint i=1; i<nic-1; i++){
		uint j = 0;
		q_chi = (q[i+skipi][j+skipj+1] - q[i+skipi][j+skipj]);
		q_eta = (q[i+skipi+1][j+skipj] - q[i+skipi-1][j+skipj])/2.0;
		grad_q[i][j][0] = chi_x[i][j]*q_chi + eta_x[i][j]*q_eta;
		grad_q[i][j][1] = chi_y[i][j]*q_chi + eta_y[i][j]*q_eta;
		j = njc-1;
		q_chi = (q[i+skipi][j+skipj] - q[i+skipi][j+skipj-1]);
		q_eta = (q[i+skipi+1][j+skipj] - q[i+skipi-1][j+skipj])/2.0;
		grad_q[i][j][0] = chi_x[i][j]*q_chi + eta_x[i][j]*q_eta;
		grad_q[i][j][1] = chi_y[i][j]*q_chi + eta_y[i][j]*q_eta;
	}

	for(uint j=1; j<njc-1; j++){
		uint i = 0;
		q_chi = (q[i+skipi][j+skipj+1] - q[i+skipi][j+skipj-1])/2.0; 
		q_eta = (q[i+skipi+1][j+skipj] - q[i+skipi][j+skipj]);
		grad_q[i][j][0] = chi_x[i][j]*q_chi + eta_x[i][j]*q_eta;
		grad_q[i][j][1] = chi_y[i][j]*q_chi + eta_y[i][j]*q_eta;
		i = nic-1;
		q_chi = (q[i+skipi][j+skipj+1] - q[i+skipi][j+skipj-1])/2.0; 
		q_eta = (q[i+skipi][j+skipj] - q[i+skipi-1][j+skipj]);
		grad_q[i][j][0] = chi_x[i][j]*q_chi + eta_x[i][j]*q_eta;
		grad_q[i][j][1] = chi_y[i][j]*q_chi + eta_y[i][j]*q_eta;
	}

	{
		uint i=0; uint j=0;
		q_chi = (q[i+skipi][j+skipj+1] - q[i+skipi][j+skipj]);
		q_eta = (q[i+skipi+1][j+skipj] - q[i+skipi][j+skipj]);
		grad_q[i][j][0] = chi_x[i][j]*q_chi + eta_x[i][j]*q_eta;
		grad_q[i][j][1] = chi_y[i][j]*q_chi + eta_y[i][j]*q_eta;
		i=nic-1; j=0;
		q_chi = (q[i+skipi][j+skipj+1] - q[i+skipi][j+skipj]);
		q_eta = (q[i+skipi][j+skipj] - q[i+skipi-1][j+skipj]);
		grad_q[i][j][0] = chi_x[i][j]*q_chi + eta_x[i][j]*q_eta;
		grad_q[i][j][1] = chi_y[i][j]*q_chi + eta_y[i][j]*q_eta;
		i=nic-1; j=njc-1;
		q_chi = (q[i+skipi][j+skipj] - q[i+skipi][j+skipj-1]);
		q_eta = (q[i+skipi][j+skipj] - q[i+skipi-1][j+skipj]);
		grad_q[i][j][0] = chi_x[i][j]*q_chi + eta_x[i][j]*q_eta;
		grad_q[i][j][1] = chi_y[i][j]*q_chi + eta_y[i][j]*q_eta;
		i=0; j=njc-1;
		q_chi = (q[i+skipi][j+skipj] - q[i+skipi][j+skipj-1]);
		q_eta = (q[i+skipi+1][j+skipj] - q[i+skipi][j+skipj]);
		grad_q[i][j][0] = chi_x[i][j]*q_chi + eta_x[i][j]*q_eta;
		grad_q[i][j][1] = chi_y[i][j]*q_chi + eta_y[i][j]*q_eta;
	}
}

template<class Tx>
Mesh<Tx>::Mesh(std::shared_ptr<Config<Tx>> val_config){
	config = val_config;
	ni = config->geometry->ni;
	nj = config->geometry->nj;
	j1 = config->geometry->tail;
	nb = ni - 2*j1 + 1;
	nic = ni - 1;
	njc = nj - 1;
	
	xv = Array2D<Tx>(ni, nj);
	yv = Array2D<Tx>(ni, nj);

	if(config->geometry->format == "simple"){
		simple_loader(config->geometry->filename);
	}
	else if(config->geometry->format == "p3d"){
		plot3d_loader(config->geometry->filename);
	}
	else{
		spdlog::get("console")->critical("file format not found!");
	}
	xc = Array2D<Tx>(nic, njc);
	yc = Array2D<Tx>(nic, njc);

	volume = Array2D<Tx>(nic, njc);

	normal_eta = Array3D<Tx>(ni-1, nj, 2U);
	normal_chi = Array3D<Tx>(ni, nj-1, 2U);

	ds_eta = Array2D<Tx>(nic, njc);
	ds_chi = Array2D<Tx>(nic, njc);

	chi_x = Array2D<Tx>(nic, njc);
	eta_x = Array2D<Tx>(nic, njc);

	chi_y = Array2D<Tx>(nic, njc);
	eta_y = Array2D<Tx>(nic, njc);

	x_chi = Array2D<Tx>(nic, njc);
	x_eta = Array2D<Tx>(nic, njc);

	y_chi = Array2D<Tx>(nic, njc);
	y_eta = Array2D<Tx>(nic, njc);
	
};

template<class Tx>
void Mesh<Tx>::setup(){
	calc_metrics();
	fluid_model = std::make_shared<FluidModel<Tx, adouble>>(config->freestream->p_inf, config->freestream->rho_inf,
															config->freestream->T_inf, config->freestream->mu_inf,
															config->freestream->pr_inf);
	solution = std::make_shared<Solution<Tx>>(this->shared_from_this());
	equation = std::make_shared<EulerEquation<Tx, adouble>>(this->shared_from_this(), config);
	iomanager = std::make_shared<IOManager<Tx>>(this->shared_from_this(), config);
	equation->initialize();

#if defined(ENABLE_ARMA)
	linearsolver = std::make_shared<LinearSolverArma<Tx>>(this->shared_from_this(), config);
#endif
	
#if defined(ENABLE_EIGEN)
	linearsolver = std::make_shared<LinearSolverEigen<Tx>>(this->shared_from_this(), config);
#endif
	
#if defined(ENABLE_PETSC)
	linearsolver = std::make_shared<LinearSolverPetsc<Tx>>(this->shared_from_this(), config);
#endif

}

template<class Tx>
Mesh<Tx>::Mesh(std::shared_ptr<Mesh<Tx>> mesh, const uint nskipi, const uint nskipj, const uint refine){
	config = mesh->config;
	if(!refine){
		ni = std::ceil((float)mesh->ni/(nskipi+1U));
		nj = std::ceil((float)mesh->nj/(nskipj+1U));
		j1 = std::ceil((float)mesh->j1/(nskipi+1U));
		nb = std::ceil((float)mesh->nb/(nskipi+1U));
	} else{
		ni = mesh->ni + (mesh->ni-1)*nskipi;
		nj = mesh->nj + (mesh->nj-1)*nskipj;
		j1 = mesh->j1 + (mesh->j1-1)*nskipi;
		nb = mesh->nb + (mesh->nb-1)*nskipi + 1;
		std::cout<<"nb "<<nb<<std::endl;
		std::cout<<"j1 "<<j1<<std::endl;
		//spdlog::get("console")->debug("NB = {}.", nb);
		//spdlog::get("console")->debug("j1 = {}.", j1);
	}
	
	nic = ni - 1;
	njc = nj - 1;
	
	xv = Array2D<Tx>(ni, nj);
	yv = Array2D<Tx>(ni, nj);

	
	for(uint i=0; i<ni; i++){
		for(uint j=0; j<nj; j++){
			if(!refine){
				uint io = i*(nskipi + 1U);
				uint jo = j*(nskipj + 1U);
				xv[i][j] = mesh->xv[io][jo];
				yv[i][j] = mesh->yv[io][jo];
			}
		}
	}

	if(refine){
		for(uint i=0; i<ni; i+=2){
			for(uint j=0; j<nj; j+=2){
				xv[i][j] = mesh->xv[i/2][j/2];
				yv[i][j] = mesh->yv[i/2][j/2];
			}
		}

		for(uint i=1; i<ni; i+=2){
			for(uint j=0; j<nj; j+=1){
				xv[i][j] = 0.5*(xv[i+1][j] + xv[i-1][j]);
				yv[i][j] = 0.5*(yv[i+1][j] + yv[i-1][j]);
			}
		}
		
		for(uint i=0; i<ni; i+=1){
			for(uint j=1; j<nj; j+=2){
				xv[i][j] = 0.5*(xv[i][j+1] + xv[i][j-1]);
				yv[i][j] = 0.5*(yv[i][j+1] + yv[i][j-1]);
			}
		}

	}
	xc = Array2D<Tx>(nic, njc);
	yc = Array2D<Tx>(nic, njc);

	volume = Array2D<Tx>(nic, njc);

	normal_eta = Array3D<Tx>(ni-1, nj, 2U);
	normal_chi = Array3D<Tx>(ni, nj-1, 2U);

	ds_eta = Array2D<Tx>(nic, njc);
	ds_chi = Array2D<Tx>(nic, njc);

	chi_x = Array2D<Tx>(nic, njc);
	eta_x = Array2D<Tx>(nic, njc);

	chi_y = Array2D<Tx>(nic, njc);
	eta_y = Array2D<Tx>(nic, njc);

	x_chi = Array2D<Tx>(nic, njc);
	x_eta = Array2D<Tx>(nic, njc);

	y_chi = Array2D<Tx>(nic, njc);
	y_eta = Array2D<Tx>(nic, njc);
	
	solution = std::make_shared<Solution<Tx>>(this);
	calc_metrics();
};

template<class Tx>
Mesh<Tx>::~Mesh(){
	
}


#endif
