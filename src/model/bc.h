#ifndef _BC_H
#define _BC_H
#include "common.h"

#define bottom 0
#define right 1
#define top 2
#define left 3

template<class T, class Tad>
class BoundaryCondition{
public:
	virtual void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p){};
};

template<class T, class Tad>
class BoundaryConditionFreestream: public BoundaryCondition<T, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	T rho_inf, u_inf, v_inf, p_inf;
	BoundaryConditionFreestream(std::string val_name, std::shared_ptr<Mesh<T>> mesh, std::shared_ptr<Config<T>> config,
								uint val_face, uint val_start, uint val_end): BoundaryCondition<T, Tad>(){
		name = val_name;
		face = val_face;
		start = val_start;
		end = val_end;
		rho_inf = config->freestream->rho_inf;
		u_inf = config->freestream->u_inf;
		v_inf = config->freestream->v_inf;
		p_inf = config->freestream->p_inf;
		ni = mesh->ni;
		nj = mesh->nj;
		nic = ni - 1;
		njc = nj - 1;
	};
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p){
		uint iend, jend;
		if(face == bottom || face == top){
			if(face == bottom)
				jend = 0;
			else
				jend = njc + 1;

#pragma omp parallel for
			for(uint i=start; i<=end; i++){
				rho[i][jend] = rho_inf;
				u[i][jend] = u_inf;
				v[i][jend] = v_inf;
				p[i][jend] = p_inf;
			}
		}

		if(face == left || face == right){
			if(face == left)
				iend = 0;
			else
				iend = nic + 1;

#pragma omp parallel for
			for(uint j=start; j<=end; j++){
				rho[iend][j] = rho_inf;
				u[iend][j] = u_inf;
				v[iend][j] = v_inf;
				p[iend][j] = p_inf;
			}
		}
		
	};
};

template<class T, class Tad>
class BoundaryConditionInviscidWall: public BoundaryCondition<T, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	T rho_inf, u_inf, v_inf, p_inf;
	std::shared_ptr<Mesh<T>> mesh;
	BoundaryConditionInviscidWall(std::string val_name, std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config<T>> config,
								  uint val_face, uint val_start, uint val_end): BoundaryCondition<T, Tad>(){
		name = val_name;
		mesh = val_mesh;
		face = val_face;
		start = val_start;
		end = val_end;
		rho_inf = config->freestream->rho_inf;
		u_inf = config->freestream->u_inf;
		v_inf = config->freestream->v_inf;
		p_inf = config->freestream->p_inf;
		ni = mesh->ni;
		nj = mesh->nj;
		nic = ni - 1;
		njc = nj - 1;
	};
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p){
		uint iend, jend;
		if(face == bottom || face == top){
			if(face == bottom)
				jend = 0;
			else
				jend = njc + 1;

#pragma omp parallel for
			for(uint i=start; i<=end; i++){
				Tad un, ds;
				ds = mesh->normal_eta[i-1][jend][0]*mesh->normal_eta[i-1][jend][0] +
					mesh->normal_eta[i-1][jend][1]*mesh->normal_eta[i-1][jend][1];
				
				p[i][jend] = 1.5*p[i][1] - 0.5*p[i][2];
				rho[i][jend] = 1.5*rho[i][1] - 0.5*rho[i][2];
				un = u[i][1]*mesh->normal_eta[i-1][jend][0] + v[i][1]*mesh->normal_eta[i-1][jend][1];
				u[i][jend] = u[i][1] - 2.0*un*mesh->normal_eta[i-1][jend][0]/ds;
				v[i][jend] = v[i][1] - 2.0*un*mesh->normal_eta[i-1][jend][1]/ds;
			}
		}

		if(face == left || face == right){
			if(face == left)
				iend = 0;
			else
				iend = nic + 1;

#pragma omp parallel for
			for(uint j=start; j<=end; j++){
				spdlog::get("console")->critical("Boundary condition not implemented!");
			}
		}
		
	};
};

template<class T, class Tad>
class BoundaryConditionViscousWall: public BoundaryCondition<T, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	T rho_inf, u_inf, v_inf, p_inf;
	std::shared_ptr<Mesh<T>> mesh;
	BoundaryConditionViscousWall(std::string val_name, std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config<T>> config,
								  uint val_face, uint val_start, uint val_end): BoundaryCondition<T, Tad>(){
		name = val_name;
		mesh = val_mesh;
		face = val_face;
		start = val_start;
		end = val_end;
		rho_inf = config->freestream->rho_inf;
		u_inf = config->freestream->u_inf;
		v_inf = config->freestream->v_inf;
		p_inf = config->freestream->p_inf;
		ni = mesh->ni;
		nj = mesh->nj;
		nic = ni - 1;
		njc = nj - 1;
	};
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p){
		uint iend, jend;
		if(face == bottom || face == top){
			if(face == bottom)
				jend = 0;
			else
				jend = njc + 1;

#pragma omp parallel for
			for(uint i=start; i<=end; i++){
				Tad un, ds;
				ds = mesh->normal_eta[i-1][jend][0]*mesh->normal_eta[i-1][jend][0] +
					mesh->normal_eta[i-1][jend][1]*mesh->normal_eta[i-1][jend][1];
				if(face == bottom){
					p[i][jend] = 1.5*p[i][1] - 0.5*p[i][2];
					rho[i][jend] = 1.5*rho[i][1] - 0.5*rho[i][2];
					u[i][jend] = -u[i][1];
					v[i][jend] = -v[i][1];
				}
				else{
					//spdlog::get("console")->debug("{} {} {}", i, start, end);
					p[i][jend] = 1.5*p[i][njc] - 0.5*p[i][njc-1];
					rho[i][jend] = 1.5*rho[i][njc] - 0.5*rho[i][njc-1];
					u[i][jend] = -u[i][njc];
					v[i][jend] = -v[i][njc];
				}
			}
		}

		if(face == left || face == right){
			if(face == left)
				iend = 0;
			else
				iend = nic + 1;

#pragma omp parallel for
			for(uint j=start; j<=end; j++){
				spdlog::get("console")->critical("Boundary condition not implemented!");
			}
		}
		
	};
};


template<class T, class Tad>
class BoundaryConditionWake: public BoundaryCondition<T, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	T rho_inf, u_inf, v_inf, p_inf;
	std::shared_ptr<Mesh<T>> mesh;
	BoundaryConditionWake(std::string val_name, std::shared_ptr<Mesh<T>> val_mesh, std::shared_ptr<Config<T>> config,
						  uint val_face, uint val_start, uint val_end): BoundaryCondition<T, Tad>(){
		name = val_name;
		mesh = val_mesh;
		face = val_face;
		start = val_start;
		end = val_end;
		rho_inf = config->freestream->rho_inf;
		u_inf = config->freestream->u_inf;
		v_inf = config->freestream->v_inf;
		p_inf = config->freestream->p_inf;
		ni = mesh->ni;
		nj = mesh->nj;
		nic = ni - 1;
		njc = nj - 1;
	};
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p){
		uint iend, jend;
		if(face == bottom || face == top){
			if(face == bottom)
				jend = 0;
			else
				jend = njc + 1;

#pragma omp parallel for
			for(uint i=start; i<=end; i++){
				rho[i][jend] = rho[nic+1-i][1];
				u[i][jend] = u[nic+1-i][1];
				v[i][jend] = v[nic+1-i][1];
				p[i][jend] = p[nic+1-i][1];
			}	

#pragma omp parallel for
			for(uint i=start; i<=end; i++){
				rho[nic+1-i][jend] = rho[i][1];
				u[nic+1-i][jend] = u[i][1];
				v[nic+1-i][jend] = v[i][1];
				p[nic+1-i][jend] = p[i][1];
			}
		}

		if(face == left || face == right){
			if(face == left)
				iend = 0;
			else
				iend = nic + 1;

#pragma omp parallel for
			for(uint j=start; j<=end; j++){
				spdlog::get("console")->critical("Boundary condition not implemented!");
			}
		}
		
	};
};
#endif
