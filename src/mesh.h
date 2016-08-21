#ifndef _MESH_H
#define _MESH_H

#include "common.h"
#include "utils.h"
#include "flux.h"
#include "solution.h"
#include "config.h"

template<class T>
class Solution;

template <class T>
class Mesh{
public:
	Config *config;
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
	Mesh(Config *config);
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
Mesh<T>::Mesh(Config *val_config){
	config = val_config;
	ni = config->geometry->ni;
	nj = config->geometry->nj;
	j1 = config->geometry->tail;
	nb = ni - 2*j1 + 1;
	nic = ni - 1;
	njc = nj - 1;
	
	xv = allocate_2d_array<T>(ni, nj);
	yv = allocate_2d_array<T>(ni, nj);

	simple_loader(config->geometry->filename);
	
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
	
	xv = allocate_2d_array<T>(ni, nj);
	yv = allocate_2d_array<T>(ni, nj);

	
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

	const unsigned int shapetmp[] = {5U, 2U};
	double *xc_array = new double[nic*njc];
	double *yc_array = new double[nic*njc];
	double *q_array = new double[nic*njc*4];

	
	for(int i=0; i<nic; i++){
		for(int j=0; j<njc; j++){
			xc_array[i*njc + j] = xc[i][j];
			yc_array[i*njc + j] = yc[i][j];
			for(int k=0; k<4; k++){
				q_array[i*njc*4 + j*4 + k] = q[i][j][k];
			}
		}
	}
	cnpy::npz_save(filename,"xc",xc_array,shape,2,"w");
	cnpy::npz_save(filename,"yc",yc_array,shape,2,"a");
	cnpy::npz_save(filename,"q",q_array,shapeq,3,"a");


	delete[] xc_array;
	delete[] yc_array;
	delete[] q_array;
}


template <class T>
void write_restart_file(const Mesh<T> *mesh, std::string filename){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	T ***q = mesh->solution->q;
	T **xc = mesh->xc;
	T **yc = mesh->yc;
	std::ofstream outfile(filename,std::ofstream::binary);
	for(int i=0; i<nic; i++){
		for(int j=0; j<njc; j++){
			outfile.write(reinterpret_cast<const char*>(q[i][j]), sizeof(T)*nq);
		}
	}
	outfile.close();
	
}

template <class T>
void read_restart_file(const Mesh<T> *mesh, std::string filename){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	T ***q = mesh->solution->q;
	T **xc = mesh->xc;
	T **yc = mesh->yc;
	std::ifstream infile(filename,std::ofstream::binary);

	infile.seekg (0,infile.end);
	long size = infile.tellg();
	long size_expected = nic*njc*nq*sizeof(T);
	infile.seekg (0);
	assert(size == size_expected);
	
	for(int i=0; i<nic; i++){
		for(int j=0; j<njc; j++){
			infile.read(reinterpret_cast<char*>(q[i][j]), sizeof(T)*nq);
		}
	}
	infile.close();
	
}

template <class T>
void write_surface_file(const Mesh<T> *mesh, std::string filename){
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	T ***q = mesh->solution->q;
	T **xc = mesh->xc;
	T **yc = mesh->yc;
	double p_inf = 1/1.4;
	double rho_inf = 1.0;
	double u_inf = 0.5;
	int j1 = mesh->j1-1;
	double rho, u, v, p, x, cp;
	std::ofstream outfile;
	outfile.open(filename);

	for(uint i=j1; i<j1+mesh->nb; i++){
		primvars<double>(q[i][0], &rho, &u, &v, &p);
		x = xc[i][0];
		cp = (p - p_inf)/(0.5*rho_inf*u_inf*u_inf);
		outfile<<x<<" "<<cp<<std::endl;
	}
	outfile.close();
}

#endif
