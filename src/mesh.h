#ifndef _MESH_H
#define _MESH_H

#include "common.h"
#include "utils.h"
#include "flux.h"
#include "solution.h"
template<class T>
class Solution;

template <class T>
class Mesh{
public:
	uint ni, nj;
	uint nic, njc;
	uint j1, nb;
	T **xv, **yv;
	T **xc, **yc;

	T ***normal_eta;
	T ***normal_chi;

	T **ds_eta;
	T **ds_chi;

	T **volume;

	Solution<T> *solution;
	
public:
	Mesh();
	Mesh(const Mesh<T>* mesh, const uint nskipi=0, const uint nskipj=0, const uint refine=0);
	~Mesh();
	void calc_metrics();
	void simple_loader(std::string filename);
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
void Mesh<T>::calc_metrics(){
	T ***reta = allocate_3d_array<T>(nic, nj, 2U);
	T ***rchi = allocate_3d_array<T>(ni, njc, 2U);

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

	release_3d_array<T>(reta, nic, nj, 2U);
	release_3d_array<T>(rchi, ni, njc, 2U);
	
};

template<class T>
Mesh<T>::Mesh(){
	ni = 281;
	nj = 51;
	j1 = 41;
	nb = 200;
	nic = ni - 1;
	njc = nj - 1;
	
	xv = allocate_2d_array<T>(ni, nj);
	yv = allocate_2d_array<T>(ni, nj);

	simple_loader("grid.unf2");
	
	xc = allocate_2d_array<T>(nic, njc);
	yc = allocate_2d_array<T>(nic, njc);

	volume = allocate_2d_array<T>(nic, njc);

	normal_eta = allocate_3d_array<T>(ni-1, nj, 2U);
	normal_chi = allocate_3d_array<T>(ni, nj-1, 2U);

	ds_eta = allocate_2d_array<T>(nic, njc);
	ds_chi = allocate_2d_array<T>(nic, njc);

	solution = new Solution<T>(this);
	this->calc_metrics();
};


template<class T>
Mesh<T>::Mesh(const Mesh<T>* mesh, const uint nskipi, const uint nskipj, const uint refine){
	if(!refine){
		ni = std::ceil((float)mesh->ni/(nskipi+1U));
		nj = std::ceil((float)mesh->nj/(nskipj+1U));
		j1 = std::ceil((float)mesh->j1/(nskipi+1U));
		nb = std::ceil((float)mesh->nb/(nskipi+1U));
	} else{
		ni = mesh->ni + (mesh->ni-1)*nskipi;
		nj = mesh->nj + (mesh->nj-1)*nskipj;
		j1 = mesh->j1 + (mesh->j1-1)*nskipi;
		nb = mesh->nb + (mesh->nb-1)*nskipi;
	}
	
	nic = ni - 1;
	njc = nj - 1;
	
	xv = allocate_2d_array<T>(ni, nj);
	yv = allocate_2d_array<T>(ni, nj);

	
	for(uint i=0; i<ni; i++){
		for(uint j=0; j<nj; j++){
			if(!refine){
				uint io = i*(nskipi + 1U);
				uint jo = j*(nskipj + 1U);
				xv[i][j] = mesh->xv[io][jo];
				yv[i][j] = mesh->yv[io][jo];
			} else{
				uint i1 = std::floor((float)i/(nskipi + 1U));
				uint j1 = std::floor((float)j/(nskipj + 1U));
				uint i2 = std::ceil((float)i/(nskipi + 1U));
				uint j2 = std::ceil((float)j/(nskipj + 1U));
				xv[i][j] = 0.5*(mesh->xv[i1][j1] + mesh->xv[i2][j2]);
				yv[i][j] = 0.5*(mesh->yv[i1][j1] + mesh->yv[i2][j2]);
			}
		}
	}
	
	xc = allocate_2d_array<T>(nic, njc);
	yc = allocate_2d_array<T>(nic, njc);

	volume = allocate_2d_array<T>(nic, njc);

	normal_eta = allocate_3d_array<T>(ni-1, nj, 2U);
	normal_chi = allocate_3d_array<T>(ni, nj-1, 2U);

	ds_eta = allocate_2d_array<T>(nic, njc);
	ds_chi = allocate_2d_array<T>(nic, njc);

	solution = new Solution<T>(this);
	//solution = new Solution<T>(this, mesh, nskipi, nskipj, refine);
	this->calc_metrics();
};

template<class T>
Mesh<T>::~Mesh(){
	release_2d_array(xv, ni, nj);
	release_2d_array(yv, ni, nj);

	release_2d_array(xc, nic, njc);
	release_2d_array(yc, nic, njc);
	release_2d_array(volume, nic, njc);
	release_2d_array(ds_eta, nic, njc);
	release_2d_array(ds_chi, nic, njc);
		
	
	release_3d_array(normal_eta, ni-1, nj, 2U);
	release_3d_array(normal_chi, ni, nj-1, 2U);
	
	delete solution;
}



template <class T>
void write_solution(const Mesh<T> *mesh, std::string filename){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	T ***q = mesh->solution->q;
	T **xc = mesh->xc;
	T **yc = mesh->yc;
	
	std::ofstream outfile;
	outfile.open(filename);
	char buffer [500];
	outfile<<"title = \"Solution\""<<"\n";
	outfile<<"variables = \"x\" \"y\" \"rho\" \"u\" \"v\" \"p\""<<"\n"; 
	outfile<<"zone i="<<nic<<", j="<<njc<<", f=point\n";
	for(uint j=0; j<njc; j++){
		for(uint i=0; i<nic; i++){
			sprintf(buffer, "%8.14E %8.14E %8.14E %8.14E %8.14E %8.14E\n", xc[i][j], yc[i][j], q[i][j][0], q[i][j][1], q[i][j][2], q[i][j][3]);
			outfile<<buffer;
		}
	}
	
	outfile.close();
	spdlog::get("console")->info("Writing tecplot data file {}.", filename);
}


template <class T>
void write_solution_npy(const Mesh<T> *mesh, std::string filename){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	T ***q = mesh->solution->q;
	T **xc = mesh->xc;
	T **yc = mesh->yc;
	const unsigned int shape[] = {nic, njc};
	const unsigned int shapeq[] = {nic, njc, 4};

	cnpy::npz_save(filename,"xc",xc,shape,1,"w");
	cnpy::npz_save(filename,"yc",yc,shape,1,"a");
	cnpy::npz_save(filename,"q",q,shapeq,1,"a");
	
}


#endif
