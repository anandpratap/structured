#ifndef _SOLVER_H
#define _SOLVER_H
#include "common.h"
#include "utils.h"
#include "mesh.h"

template<class T>
class Solver{
public:
	int nnz;
	unsigned int *rind = nullptr;
	unsigned int *cind = nullptr;
	double *values = nullptr;
	int options[4] = {0,0,0,0};
	
	Mesh<T> *mesh;
	void calc_residual();
	void solve();
	adouble ***a_q, ***a_rhs;
	
	Solver(Mesh<T> *val_mesh);
	~Solver();
};

template <class T>
Solver<T>::Solver(Mesh<T> *val_mesh){
	mesh = val_mesh;
	a_q = allocate_3d_array<adouble>(mesh->nic, mesh->njc, mesh->solution->nq);
	a_rhs = allocate_3d_array<adouble>(mesh->nic, mesh->njc, mesh->solution->nq);
}
template <class T>
Solver<T>::~Solver(){
	release_3d_array(a_q, mesh->nic, mesh->njc, mesh->solution->nq);
	release_3d_array(a_rhs, mesh->nic, mesh->njc, mesh->solution->nq);
}


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
	double dt = 1e2;
	double t = 0.0;
	while(1){
		counter += 1;
		trace_on(1);
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq; k++){
					mesh->solution->a_q[i*njc*nq + j*nq + k] <<= mesh->solution->q[i][j][k];
					mesh->solution->tmp_q[i*njc*nq + j*nq + k] = mesh->solution->q[i][j][k];
				}
			}
		}
		
		calc_residual();

		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq; k++){
					mesh->solution->a_rhs[i*njc*nq + j*nq + k] >>= mesh->solution->rhs[i][j][k];
				}
			}
		}
		
		trace_off();
		sparse_jac(1,nic*njc*nq,nic*njc*nq,0,mesh->solution->tmp_q,&nnz,&rind,&cind,&values,options);
		std::cout<<"NNZ = "<<nnz<<std::endl;
		for(int i=0; i<nnz; i++){
			mesh->solution->arma_jac(rind[i],cind[i]) = -values[i];
		}
		
		for(int i=0; i<nic*njc*nq; i++){
			mesh->solution->arma_jac(i,i) += 1.0/dt;
		}

		for(int i=0; i<nic*njc*nq; i++){
			mesh->solution->arma_rhs(i, 0) = mesh->solution->a_rhs[i].value();
		}

		mesh->solution->arma_dq = arma::spsolve(mesh->solution->arma_jac, mesh->solution->arma_rhs, "lapack");

		free(rind); rind=nullptr;
		free(cind); cind=nullptr;
		free(values); values=nullptr;
		
		
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				for(uint k=0; k<nq; k++){
					//q[i][j][k] = q[i][j][k] + rhs[i][j][k]*dt;
					q[i][j][k] = q[i][j][k] + mesh->solution->arma_dq(i*njc*nq + j*nq + k, 0);
				}
			}
		}
		t += dt;
		counter += 1;
		
		if (l2norm < 1e-8) break;
		if(counter % 1 == 0){
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
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			for(uint k=0; k<nq; k++){
				a_q[i][j][k] = mesh->solution->a_q[i*njc*nq + j*nq + k];
			}
		}
	}
	
	adouble ***q = a_q;
	static adouble **rho = allocate_2d_array<adouble>(nic+2, njc+2);
	static adouble **u = allocate_2d_array<adouble>(nic+2, njc+2);
	static adouble **v = allocate_2d_array<adouble>(nic+2, njc+2);
	static adouble **p = allocate_2d_array<adouble>(nic+2, njc+2);
	
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			primvars<adouble>(q[i][j], &rho[i+1][j+1], &u[i+1][j+1], &v[i+1][j+1], &p[i+1][j+1]);
			// std::cout<<rho[i][j]<<std::endl;
			// std::cout<<u[i][j]<<std::endl;
			// std::cout<<v[i][j]<<std::endl;
			// std::cout<<p[i][j]<<std::endl;
			for(uint k=0; k<nq; k++){
				a_rhs[i][j][k] = 0.0;
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


	uint j1 = mesh->j1;
	uint nb = mesh->nb;

	adouble un, ds;

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


	
	static adouble** rholft_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** ulft_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** vlft_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** plft_xi = allocate_2d_array<adouble>(ni, njc);

	static adouble** rhorht_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** urht_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** vrht_xi = allocate_2d_array<adouble>(ni, njc);
	static adouble** prht_xi = allocate_2d_array<adouble>(ni, njc);

	first_order_xi(ni, nj, rho, rholft_xi, rhorht_xi);
	first_order_xi(ni, nj, u, ulft_xi, urht_xi);
	first_order_xi(ni, nj, v, vlft_xi, vrht_xi);
	first_order_xi(ni, nj, p, plft_xi, prht_xi);


	static adouble** rholft_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** ulft_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** vlft_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** plft_eta = allocate_2d_array<adouble>(nic, nj);

	static adouble** rhorht_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** urht_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** vrht_eta = allocate_2d_array<adouble>(nic, nj);
	static adouble** prht_eta = allocate_2d_array<adouble>(nic, nj);

	first_order_eta(ni, nj, rho, rholft_eta, rhorht_eta);
	first_order_eta(ni, nj, u, ulft_eta, urht_eta);
	first_order_eta(ni, nj, v, vlft_eta, vrht_eta);
	first_order_eta(ni, nj, p, plft_eta, prht_eta);


	static adouble*** flux_xi = allocate_3d_array<adouble>(ni, njc, 4U);
	static adouble*** flux_eta = allocate_3d_array<adouble>(nic, nj, 4U);
	
	for(uint i=0; i< ni; i++){
		for(uint j=0; j< njc; j++){
			roeflux<adouble>(mesh->normal_chi[i][j][0], mesh->normal_chi[i][j][1],
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
			roeflux<adouble>(mesh->normal_eta[i][j][0], mesh->normal_eta[i][j][1],
					   rholft_eta[i][j], ulft_eta[i][j], vlft_eta[i][j], plft_eta[i][j],
					   rhorht_eta[i][j], urht_eta[i][j], vrht_eta[i][j], prht_eta[i][j],
					   flux_eta[i][j]);
		}
	}

	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] -= (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
				a_rhs[i][j][k] -= (flux_xi[i+1][j][k] - flux_xi[i][j][k]);
			}
		}
	}

	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] /= mesh->volume[i][j];
			}
		}
	}

	
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			for(uint k=0; k<nq; k++){
				mesh->solution->a_rhs[i*njc*nq + j*nq + k] = a_rhs[i][j][k];
			}
		}
	}
	
	
};

#endif
