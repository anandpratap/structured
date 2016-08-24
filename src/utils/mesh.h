#ifndef _MESH_H
#define _MESH_H
#include "common.h"
#include "utils.h"
#include "solution.h"
#include "config.h"

template<class T>
class Solution;

template <class T>
class Mesh{
 public:
	std::shared_ptr<Config<T>> config;
	uint ni, nj;
	uint nic, njc;
	uint j1, nb;
	Array2D<T> xv, yv;
	Array2D<T> xc, yc;

	Array3D<T> normal_eta;
	Array3D<T> normal_chi;

	Array2D<T> ds_eta;
	Array2D<T> ds_chi;

	Array2D<T> volume;

	std::shared_ptr<Solution<T>> solution;
	
public:
	Mesh(std::shared_ptr<Config<T>> config);
	Mesh(std::shared_ptr<Mesh<T>> mesh, const uint nskipi=0, const uint nskipj=0, const uint refine=0);
	~Mesh();
	void calc_metrics();
	void simple_loader(std::string filename);
	void plot3d_loader(std::string filename);
};

template<class T>
void Mesh<T>::simple_loader(std::string filename){
	std::ifstream infile(filename);
	infile >> std::fixed >> std::setprecision(20);
	for(uint j=0; j<nj; j++){
		for(uint i=0; i<ni; i++){  
			infile >> xv[i][j] >> yv[i][j];
		}
	}
	infile.close();
}

template<class T>
void Mesh<T>::plot3d_loader(std::string filename){
	std::ifstream infile(filename);
	infile >> std::fixed >> std::setprecision(20);
	int nblock, ni_, nj_;
	infile >> nblock;
	infile >> ni_ >> nj_;
	assert(nblock == 1);
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

template<class T>
void Mesh<T>::calc_metrics(){
	auto reta = Array3D<T>(nic, nj, 2U);
	auto rchi = Array3D<T>(ni, njc, 2U);

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

};

template<class T>
Mesh<T>::Mesh(std::shared_ptr<Config<T>> val_config){
	config = val_config;
	ni = config->geometry->ni;
	nj = config->geometry->nj;
	j1 = config->geometry->tail;
	nb = ni - 2*j1 + 1;
	nic = ni - 1;
	njc = nj - 1;
	
	xv = Array2D<T>(ni, nj);
	yv = Array2D<T>(ni, nj);

	if(config->geometry->format == "simple"){
		simple_loader(config->geometry->filename);
	}
	else if(config->geometry->format == "p3d"){
		plot3d_loader(config->geometry->filename);
	}
	else{
		spdlog::get("console")->critical("file format not found!");
	}
	xc = Array2D<T>(nic, njc);
	yc = Array2D<T>(nic, njc);

	volume = Array2D<T>(nic, njc);

	normal_eta = Array3D<T>(ni-1, nj, 2U);
	normal_chi = Array3D<T>(ni, nj-1, 2U);

	ds_eta = Array2D<T>(nic, njc);
	ds_chi = Array2D<T>(nic, njc);

	solution = std::make_shared<Solution<T>>(this);
	calc_metrics();
};


template<class T>
Mesh<T>::Mesh(std::shared_ptr<Mesh<T>> mesh, const uint nskipi, const uint nskipj, const uint refine){
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
	
	xv = Array2D<T>(ni, nj);
	yv = Array2D<T>(ni, nj);

	
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
	xc = Array2D<T>(nic, njc);
	yc = Array2D<T>(nic, njc);

	volume = Array2D<T>(nic, njc);

	normal_eta = Array3D<T>(ni-1, nj, 2U);
	normal_chi = Array3D<T>(ni, nj-1, 2U);

	ds_eta = Array2D<T>(nic, njc);
	ds_chi = Array2D<T>(nic, njc);

	solution = std::make_shared<Solution<T>>(this);
	calc_metrics();
};

template<class T>
Mesh<T>::~Mesh(){
	
}


#endif
