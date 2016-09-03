#ifndef _BC_H
#define _BC_H
#include "common.h"

#define bottom 0
#define right 1
#define top 2
#define left 3


template<class Tx, class Tad>
class BoundaryCondition{
public:
	virtual void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){};
};

template<class Tx, class Tad>
class BoundaryConditionFreestream: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	BoundaryConditionFreestream(std::string val_name, std::shared_ptr<Mesh<Tx>> mesh, std::shared_ptr<Config<Tx>> config,
								uint val_face, uint val_start, uint val_end): BoundaryCondition<Tx, Tad>(){
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
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
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

template<class Tx, class Tad>
class BoundaryConditionInviscidWall: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	std::shared_ptr<Mesh<Tx>> mesh;
	BoundaryConditionInviscidWall(std::string val_name, std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> config,
								  uint val_face, uint val_start, uint val_end): BoundaryCondition<Tx, Tad>(){
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
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
		uint iend, jend;
		if(face == bottom || face == top){
			if(face == bottom)
				jend = 0;
			else
				jend = njc + 1;

#pragma omp parallel for
			for(uint i=start; i<=end; i++){
				Tad un, ds;
				if(jend == 0) {
					auto nx =  mesh->normal_eta[i-1][jend][0];
					auto ny =  mesh->normal_eta[i-1][jend][1];

					ds = nx*nx + ny*ny;
					p[i][jend] = 1.5*p[i][1] - 0.5*p[i][2];
					rho[i][jend] = 1.5*rho[i][1] - 0.5*rho[i][2];
					un = u[i][1]*nx + v[i][1]*ny;
					u[i][jend] = u[i][1] - 2.0*un*nx/ds;
					v[i][jend] = v[i][1] - 2.0*un*ny/ds;

				}
				else{
					auto nx =  mesh->normal_eta[i-1][jend-1][0];
					auto ny =  mesh->normal_eta[i-1][jend-1][1];
					
					ds = nx*nx + ny*ny;
					p[i][jend] = 1.5*p[i][njc] - 0.5*p[i][njc-1];
					rho[i][jend] = 1.5*rho[i][njc] - 0.5*rho[i][njc-1];
					un = u[i][njc]*nx + v[i][njc]*ny;
					u[i][jend] = u[i][njc] - 2.0*un*nx/ds;
					v[i][jend] = v[i][njc] - 2.0*un*ny/ds;
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

template<class Tx, class Tad>
class BoundaryConditionViscousWall: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	double u_bc, v_bc;
	std::shared_ptr<Mesh<Tx>> mesh;
	BoundaryConditionViscousWall(std::string val_name, std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> config,
								 uint val_face, uint val_start, uint val_end, double val_u_bc, double val_v_bc): BoundaryCondition<Tx, Tad>(){
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
		
		u_bc = val_u_bc;
		v_bc = val_v_bc;
	};
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
		uint iend, jend;
		if(face == bottom || face == top){
			if(face == bottom)
				jend = 0;
			else
				jend = njc + 1;

#pragma omp parallel for
			for(uint i=start; i<=end; i++){
				if(face == bottom){
					p[i][jend] = 1.5*p[i][1] - 0.5*p[i][2];
					rho[i][jend] = 1.5*rho[i][1] - 0.5*rho[i][2];
					u[i][jend] = u_bc - (1.5*u[i][1] - 0.5*u[i][2]);
					v[i][jend] = v_bc - (1.5*v[i][1] - 0.5*v[i][2]);
					T[i][jend] = p[i][jend]/rho[i][jend];
				}
				else{
					p[i][jend] = 1.5*p[i][njc] - 0.5*p[i][njc-1];
					rho[i][jend] = 1.5*rho[i][njc] - 0.5*rho[i][njc-1];
					u[i][jend] = u_bc - (1.5*u[i][njc] - 0.5*u[i][njc-1]);
					v[i][jend] = v_bc - (1.5*v[i][njc] - 0.5*v[i][njc-1]);
					T[i][jend] = p[i][jend]/rho[i][jend];
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


template<class Tx, class Tad>
class BoundaryConditionWake: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	std::shared_ptr<Mesh<Tx>> mesh;
	BoundaryConditionWake(std::string val_name, std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> config,
						  uint val_face, uint val_start, uint val_end): BoundaryCondition<Tx, Tad>(){
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
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
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










template<class Tx, class Tad>
class BoundaryConditionOutflow: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	BoundaryConditionOutflow(std::string val_name, std::shared_ptr<Mesh<Tx>> mesh, std::shared_ptr<Config<Tx>> config,
								uint val_face, uint val_start, uint val_end): BoundaryCondition<Tx, Tad>(){
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
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
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
				rho[iend][j] = rho[iend-1][j];
				u[iend][j] = u[iend-1][j];
				v[iend][j] = v[iend-1][j];
				p[iend][j] = p_inf;
			}
		}
		
	};
};



template<class Tx, class Tad>
class BoundaryConditionPeriodic: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	BoundaryConditionPeriodic(std::string val_name, std::shared_ptr<Mesh<Tx>> mesh, std::shared_ptr<Config<Tx>> config,
								uint val_face, uint val_start, uint val_end): BoundaryCondition<Tx, Tad>(){
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
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
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
#pragma omp parallel for
			for(uint j=start; j<=end; j++){
				rho[0][j] = rho[nic][j];
				u[0][j] = u[nic][j];
				v[0][j] = v[nic][j];
				p[0][j] = p[nic][j];
				T[0][j] = T[nic][j];
				
				rho[nic+1][j] = rho[1][j];
				u[nic+1][j] = u[1][j];
				v[nic+1][j] = v[1][j];
				p[nic+1][j] = p[1][j];
				T[nic+1][j] = T[1][j];

				//				spdlog::get("console")->debug("periodic bc");
			}
		}
		
	};
};


template<class Tx, class Tad>
class BoundaryConditionIsothermalWall: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	uint face;
	uint start, end;
	uint ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	double u_bc, v_bc, T_bc;
	std::shared_ptr<Mesh<Tx>> mesh;
 BoundaryConditionIsothermalWall(std::string val_name, std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> config,
								 uint val_face, uint val_start, uint val_end, double val_u_bc, double val_v_bc, double val_T_bc): BoundaryCondition<Tx, Tad>(){
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

		u_bc = val_u_bc;
		v_bc = val_v_bc;
		T_bc = val_T_bc;
		
	};
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
		uint iend, jend;
		if(face == bottom || face == top){
			if(face == bottom)
				jend = 0;
			else
				jend = njc + 1;

#pragma omp parallel for
			for(uint i=start; i<=end; i++){
				if(face == bottom){
					p[i][jend] = 1.5*p[i][1] - 0.5*p[i][2];
					u[i][jend] = u_bc - (1.5*u[i][1] - 0.5*u[i][2]);
					v[i][jend] = v_bc - (1.5*v[i][1] - 0.5*v[i][2]);
					T[i][jend] = T_bc;
					rho[i][jend] = p[i][jend]/T[i][jend];
				}
				else{
					p[i][jend] = 1.5*p[i][njc] - 0.5*p[i][njc-1];
					u[i][jend] = u_bc - (1.5*u[i][njc] - 0.5*u[i][njc-1]);
					v[i][jend] = v_bc - (1.5*v[i][njc] - 0.5*v[i][njc-1]);
					T[i][jend] = T_bc;
					rho[i][jend] = p[i][jend]/T[i][jend];
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















template<class Tx, class Tad>
class BoundaryContainer{
 public:
	std::vector<BoundaryCondition<Tx, Tad>*> boundary_conditions;
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
		for(auto&& bc : boundary_conditions)
			bc->apply(rho, u, v, p, T);
	};

	~BoundaryContainer(){
		for(auto&& bc : boundary_conditions)
			delete bc;
	};
	BoundaryContainer(std::string filename, std::shared_ptr<Mesh<Tx>> & mesh, std::shared_ptr<Config<Tx>> val_config){
		auto config = cpptoml::parse_file(filename);
		auto bcs = config->get_table_array("boundary");
		for (const auto& bc : *bcs){
			std::string name = bc->get_qualified_as<std::string>("name").value_or("boundary");
			std::string type = bc->get_qualified_as<std::string>("type").value_or("");
			std::string face = bc->get_qualified_as<std::string>("face").value_or("");
			int start = bc->get_qualified_as<int64_t>("start").value_or(0);
			int end = bc->get_qualified_as<int64_t>("end").value_or(0);

			uint facei;
			if(face=="left") facei = left;
			if(face=="right") facei = right;
			if(face=="bottom") facei = bottom;
			if(face=="top") facei = top;


			double u_bc = bc->get_qualified_as<double>("u").value_or(0.0);
			double v_bc = bc->get_qualified_as<double>("v").value_or(0.0);
			double T_bc = bc->get_qualified_as<double>("T").value_or(0.0);
			
			//auto c = config->config;
			//spdlog::get("console")->debug("{}", c);
			//auto tmp_inf =  c->get_qualified_as<double>("freestream.rho_inf").value_or(1.0);
			
			if(type == "freestream"){
				auto boundarycondition = new BoundaryConditionFreestream<Tx, Tad>(name, mesh, val_config, facei, start, end);
				boundary_conditions.push_back(boundarycondition);
			}
			else if(type=="outflow"){
				auto boundarycondition = new BoundaryConditionOutflow<Tx, Tad>(name, mesh, val_config, facei, start, end);
				boundary_conditions.push_back(boundarycondition);
			}
			else if(type=="periodic"){
				auto boundarycondition = new BoundaryConditionPeriodic<Tx, Tad>(name, mesh, val_config, facei, start, end);
				boundary_conditions.push_back(boundarycondition);
			}
			else if(type == "slipwall"){
				auto boundarycondition = new BoundaryConditionInviscidWall<Tx, Tad>(name, mesh, val_config, facei, start, end);
				boundary_conditions.push_back(boundarycondition);
			}
			else if(type == "wall"){
				auto boundarycondition = new BoundaryConditionViscousWall<Tx, Tad>(name, mesh, val_config, facei, start, end, u_bc, v_bc);
				boundary_conditions.push_back(boundarycondition);
			}
			else if(type == "isothermalwall"){
				auto boundarycondition = new BoundaryConditionIsothermalWall<Tx, Tad>(name, mesh, val_config, facei, start, end, u_bc, v_bc, T_bc);
				boundary_conditions.push_back(boundarycondition);
			}
			else if(type == "wake"){
				auto boundarycondition = new BoundaryConditionWake<Tx, Tad>(name, mesh, val_config, facei, start, end);
				boundary_conditions.push_back(boundarycondition);
			}
			else{
				spdlog::get("console")->info("Wrong type of BC.");
			}
		}
	}
};
#endif
