#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "utils.h"
const double GAMMA = 1.4;
#include "flux.h"
template<class T>
void first_order_xi(uint ni, uint nj, T** q, T** ql, T** qr){
	uint njm = nj-1;
	for(uint i=0; i<ni; i++){
		for(uint j=0; j<njm; j++){
			ql[i][j] = q[i][j+1];
			qr[i][j] = q[i+1][j+1];
			//std::cout<<i<<" "<<j<<" "<<ql[i][j]<<" "<<qr[i][j]<<std::endl;

		}
	}
}

template<class T>
void first_order_eta(uint ni, uint nj, T** q, T** ql, T** qr){
	uint nim = ni-1;
  	for(uint i=0; i<nim; i++){
		for(uint j=0; j<nj; j++){
			ql[i][j] = q[i+1][j];
			qr[i][j] = q[i+1][j+1];
		}
	}
}

template<class T>
void primvars(const T Q[4], T *rho, T *u, T *v, T *p){
  T tmp_rho = Q[0];
  T tmp_u = Q[1]/tmp_rho;
  T tmp_v = Q[2]/tmp_rho;
  *rho = Q[0];
  *u = Q[1]/Q[0];
  *v = Q[2]/Q[0];
  *p = (Q[3] - 0.5*tmp_rho*(tmp_u*tmp_u + tmp_v*tmp_v))*(GAMMA-1);
}


typedef unsigned int uint;

template<class T>
class Mesh;

template<class T>
class Solution{
public:
	uint nic, njc, nq, naux;
	T ***q;
	T ***rhs;
	T **lhs;
	T ***q_aux;
public:
	Solution(const Mesh<T>* mesh);
	Solution(const Mesh<T>* mesh, const Mesh<T>* old_mesh, const uint nskipi=0, const uint nskipj=0, const uint refine=0);

	~Solution();
};
template<class T>
Solution<T>::Solution(const Mesh<T>* mesh){
	nq = 4;
	naux = 10;
	nic = mesh->nic;
	njc = mesh->njc;
	q = allocate_3d_array<T>(nic, njc, nq);
	rhs = allocate_3d_array<T>(nic, njc, nq);
	q_aux = allocate_3d_array<T>(nic, njc, naux);

	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			q[i][j][0] = 1.0;
			q[i][j][1] = 0.8;
			q[i][j][2] = 0.0;
			q[i][j][3] = 2.1057142857142863;
			for(uint k=0; k<nq; k++){
				rhs[i][j][k] = 0.0;
			}
		}
	}
}




template<class T>
Solution<T>::Solution(const Mesh<T>* mesh, const Mesh<T>* old_mesh, const uint nskipi, const uint nskipj, const uint refine): Solution<T>(mesh){
	std::cout<<nq<<std::endl;
	// interpolate
}


template<class T>
Solution<T>::~Solution(){
	release_3d_array(rhs, nic, njc, nq);
	release_3d_array(q, nic, njc, nq);
	release_3d_array(q_aux, nic, njc, naux);
}

template <class T>
class Mesh{
public:
	uint ni, nj;
	uint nic, njc;
	
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
	} else{
		ni = mesh->ni + (mesh->ni-1)*nskipi;
		nj = mesh->nj + (mesh->nj-1)*nskipj;
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

	solution = new Solution<T>(this, mesh, nskipi, nskipj, refine);
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
	printf("Data dumped..\n");
}

template<class T>
class Solver{
public:
	Mesh<T> *mesh;
	void calc_residual();
	void solve();
};
template <class T>
void Solver<T>::solve(){
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;
	double l2norm = 1e10;

	uint counter = 0;
	T ***q = mesh->solution->q;
	T ***rhs = mesh->solution->rhs;
	double dt = 1e-5;
	double t = 0.0;
	while(1){
		counter += 1;
		calc_residual();
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq; k++){
					q[i][j][k] = q[i][j][k] + rhs[i][j][k]/mesh->volume[i][j]*dt;
				}
			}
		}
		t += dt;
		counter += 1;
		
		if (l2norm < 1e-8) break;
		if(counter % 1000 == 0){
			l2norm = 0.0;
			for(uint i=0; i<nic; i++){
				for(uint j=0; j<njc; j++){
					l2norm += rhs[i][j][0]*rhs[i][j][0];
				}
			}
			l2norm = sqrt(l2norm);
			std::cout<<"l2norm = "<<l2norm;
			std::cout<<" t = "<<t;
			std::cout<<" counter = "<<counter<<std::endl;
			write_solution(mesh, "base.tec");
		}		
	}
}
template <class T>
void Solver<T>::calc_residual(){
	//std::cout<<"calc_res"<<std::endl;
	uint ni = mesh->ni;
	uint nj = mesh->nj;
	uint nic = mesh->nic;
	uint njc = mesh->njc;
	uint nq = mesh->solution->nq;

	T ***q = mesh->solution->q;
	static T **rho = allocate_2d_array<T>(nic+2, njc+2);
	static T **u = allocate_2d_array<T>(nic+2, njc+2);
	static T **v = allocate_2d_array<T>(nic+2, njc+2);
	static T **p = allocate_2d_array<T>(nic+2, njc+2);
	
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			primvars<T>(q[i][j], &rho[i+1][j+1], &u[i+1][j+1], &v[i+1][j+1], &p[i+1][j+1]);
			// std::cout<<rho[i][j]<<std::endl;
			// std::cout<<u[i][j]<<std::endl;
			// std::cout<<v[i][j]<<std::endl;
			// std::cout<<p[i][j]<<std::endl;
			for(uint k=0; k<nq; k++){
				mesh->solution->rhs[i][j][k] = 0.0;
			}
		}
	}

	for(uint i=0; i<nic+2; i++){
		rho[i][njc+1] = 1.0;
		u[i][njc+1] = 0.8;
		v[i][njc+1] = 0.0;
		p[i][njc+1] = 1.0/GAMMA;
		
	}

	for(uint j=0; j<njc+2; j++){
		rho[0][j] = 1.0;
		u[0][j] = 0.8;
		v[0][j] = 0.0;
		p[0][j] = 1.0/GAMMA;
		rho[nic+1][j] = 1.0;
		u[nic+1][j] = 0.8;
		v[nic+1][j] = 0.0;
		p[nic+1][j] = 1.0/GAMMA;
	}


	uint j1 = 41;
	uint nb = 200;

	T un, ds;

	for(uint i=0; i<nb; i++){
		ds = mesh->normal_eta[j1-1+i][0][0]*mesh->normal_eta[j1-1+i][0][0] +
			mesh->normal_eta[j1-1+i][0][1]*mesh->normal_eta[j1-1+i][0][1];

		p[j1+i][0] = 1.5*p[j1+i][1] - 0.5*p[j1+i][2];
		rho[j1+i][0] = 1.5*rho[j1+i][1] - 0.5*rho[j1+i][2];
		un = u[j1+i][1]*mesh->normal_eta[j1-1+i][0][0] + v[j1+i][1]*mesh->normal_eta[j1-1+i][0][1];
		u[j1+i][0] = u[j1+i][1] - 2*un*mesh->normal_eta[j1-1+i][0][0]/ds;
		v[j1+i][0] = v[j1+i][1] - 2*un*mesh->normal_eta[j1-1+i][0][1]/ds;
	}

	for(uint i=1; i < j1; i++){
		rho[i][0] = rho[nic+1-i][1];
		u[i][0] = u[nic+1-i][1];
		v[i][0] = v[nic+1-i][1];
		p[i][0] = p[nic+1-i][1];

		rho[nic+1-i][0] = rho[i][1];
		u[nic+1-i][0] = u[i][1];
		v[nic+1-i][0] = v[i][1];
		p[nic+1-i][0] = p[i][1];
	}


	
	static T**rholft_xi = allocate_2d_array<T>(ni, njc);
	static T**ulft_xi = allocate_2d_array<T>(ni, njc);
	static T**vlft_xi = allocate_2d_array<T>(ni, njc);
	static T**plft_xi = allocate_2d_array<T>(ni, njc);

	static T**rhorht_xi = allocate_2d_array<T>(ni, njc);
	static T**urht_xi = allocate_2d_array<T>(ni, njc);
	static T**vrht_xi = allocate_2d_array<T>(ni, njc);
	static T**prht_xi = allocate_2d_array<T>(ni, njc);

	first_order_xi(ni, nj, rho, rholft_xi, rhorht_xi);
	first_order_xi(ni, nj, u, ulft_xi, urht_xi);
	first_order_xi(ni, nj, v, vlft_xi, vrht_xi);
	first_order_xi(ni, nj, p, plft_xi, prht_xi);


	static T**rholft_eta = allocate_2d_array<T>(nic, nj);
	static T**ulft_eta = allocate_2d_array<T>(nic, nj);
	static T**vlft_eta = allocate_2d_array<T>(nic, nj);
	static T**plft_eta = allocate_2d_array<T>(nic, nj);

	static T**rhorht_eta = allocate_2d_array<T>(nic, nj);
	static T**urht_eta = allocate_2d_array<T>(nic, nj);
	static T**vrht_eta = allocate_2d_array<T>(nic, nj);
	static T**prht_eta = allocate_2d_array<T>(nic, nj);

	first_order_eta(ni, nj, rho, rholft_eta, rhorht_eta);
	first_order_eta(ni, nj, u, ulft_eta, urht_eta);
	first_order_eta(ni, nj, v, vlft_eta, vrht_eta);
	first_order_eta(ni, nj, p, plft_eta, prht_eta);


	static T***flux_xi = allocate_3d_array<T>(ni, njc, 4U);
	static T***flux_eta = allocate_3d_array<T>(nic, nj, 4U);
	
	for(uint i=0; i< ni; i++){
		for(uint j=0; j< njc; j++){
			roeflux<T>(mesh->normal_chi[i][j][0], mesh->normal_chi[i][j][1],
					   rholft_xi[i][j], ulft_xi[i][j], vlft_xi[i][j], plft_xi[i][j],
					   rhorht_xi[i][j], urht_xi[i][j], vrht_xi[i][j], prht_xi[i][j],
					   flux_xi[i][j]);
			//std::cout<<flux_xi[i][j][0]<<" "<<rholft_xi[i][j]<<" "<<rhorht_xi[i][j]<<std::endl;
			//std::cout<<flux_xi[i][j][1]<<" "<<ulft_xi[i][j]<<" "<<urht_xi[i][j]<<std::endl;
			//std::cout<<flux_xi[i][j][2]<<" "<<vlft_xi[i][j]<<" "<<vrht_xi[i][j]<<std::endl;
			//std::cout<<flux_xi[i][j][3]<<" "<<plft_xi[i][j]<<" "<<prht_xi[i][j]<<std::endl;

		}
	}

	for(uint i=0; i< nic; i++){
		for(uint j=0; j< nj; j++){
			roeflux<T>(mesh->normal_eta[i][j][0], mesh->normal_eta[i][j][1],
					   rholft_eta[i][j], ulft_eta[i][j], vlft_eta[i][j], plft_eta[i][j],
					   rhorht_eta[i][j], urht_eta[i][j], vrht_eta[i][j], prht_eta[i][j],
					   flux_eta[i][j]);
		}
	}

	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				mesh->solution->rhs[i][j][k] -= (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
				mesh->solution->rhs[i][j][k] -= (flux_xi[i+1][j][k] - flux_xi[i][j][k]);
			}
		}
	}
	
};


int main(void){
	Mesh<double> m = Mesh<double>();
	std::cout<<m.ni<<std::endl;
	//Mesh<double> m = Mesh<double>();
	//std::cout<<m.ni<<std::endl;


	Solver<double> s = Solver<double>();
	s.mesh = &m;
	s.solve();
	write_solution(&m, "base.tec");
	
	return 0;
}


