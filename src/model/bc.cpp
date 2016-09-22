#ifndef _BC_CPP
#define _BC_CPP
#include "common.h"
#include "fluid.h"
#include "bc.h"

template<class Tx, class Tad>
BoundaryConditionFreestream<Tx, Tad>::BoundaryConditionFreestream(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> mesh, std::shared_ptr<Config<Tx>> config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
														 size_t val_face, size_t val_start, size_t val_end): BoundaryCondition<Tx, Tad>(){
	fluid_model = val_fluid_model;
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
	
template<class Tx, class Tad>
void BoundaryConditionFreestream<Tx, Tad>::apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
	size_t iend, jend;
	if(face == global_bottom || face == global_top){
		if(face == global_bottom)
			jend = 0;
		else
			jend = njc + 1;

#pragma omp parallel for
		for(size_t i=start; i<=end; i++){
			rho[i][jend] = rho_inf;
			u[i][jend] = u_inf;
			v[i][jend] = v_inf;
			p[i][jend] = p_inf;
			T[i][jend] = fluid_model->get_T_prho(p[i][jend], rho[i][jend]);
		}
	}

	if(face == global_left || face == global_right){
		if(face == global_left)
			iend = 0;
		else
			iend = nic + 1;

#pragma omp parallel for
		for(size_t j=start; j<=end; j++){
			rho[iend][j] = rho_inf;
			u[iend][j] = u_inf;
			v[iend][j] = v_inf;
			p[iend][j] = p_inf;
			T[iend][j] = fluid_model->get_T_prho(p[iend][j], rho[iend][j]);
		}
	}
		
};

template<class Tx, class Tad>
BoundaryConditionInviscidWall<Tx, Tad>::BoundaryConditionInviscidWall(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
															 size_t val_face, size_t val_start, size_t val_end): BoundaryCondition<Tx, Tad>(){
	fluid_model = val_fluid_model;
	name = val_name;
	mesh = val_mesh;
	face = val_face;
	start = val_start;
	end = val_end;
	ni = mesh->ni;
	nj = mesh->nj;
	nic = ni - 1;
	njc = nj - 1;
};
	
template<class Tx, class Tad>
void BoundaryConditionInviscidWall<Tx, Tad>::apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
	size_t iend, jend;
	if(face == global_bottom || face == global_top){
		if(face == global_bottom)
			jend = 0;
		else
			jend = njc + 1;

#pragma omp parallel for
		for(size_t i=start; i<=end; i++){
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
				T[i][jend] = fluid_model->get_T_prho(p[i][jend], rho[i][jend]);
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
				T[i][jend] = fluid_model->get_T_prho(p[i][jend], rho[i][jend]);
			}
		}
	}

	if(face == global_left || face == global_right){
		if(face == global_left)
			iend = 0;
		else
			iend = nic + 1;

#pragma omp parallel for
		for(size_t j=start; j<=end; j++){
			spdlog::get("console")->critical("Boundary condition not implemented!");
		}
	}
		
};


template<class Tx, class Tad>
BoundaryConditionAdiabaticWall<Tx, Tad>::BoundaryConditionAdiabaticWall(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
															   size_t val_face, size_t val_start, size_t val_end, Tx val_u_bc, Tx val_v_bc): BoundaryCondition<Tx, Tad>(){
	fluid_model = val_fluid_model;
	name = val_name;
	mesh = val_mesh;
	face = val_face;
	start = val_start;
	end = val_end;
	ni = mesh->ni;
	nj = mesh->nj;
	nic = ni - 1;
	njc = nj - 1;
		
	u_bc = val_u_bc;
	v_bc = val_v_bc;
};
	
template<class Tx, class Tad>
void BoundaryConditionAdiabaticWall<Tx, Tad>::apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
	size_t iend, jend;
	if(face == global_bottom || face == global_top){
		if(face == global_bottom)
			jend = 0;
		else
			jend = njc + 1;

#pragma omp parallel for
		for(size_t i=start; i<=end; i++){
			if(face == global_bottom){
				T[i][jend] = 1.5*T[i][1] - 0.5*T[i][2];
				rho[i][jend] = 1.5*rho[i][1] - 0.5*rho[i][2];
				u[i][jend] = 2.0*u_bc - (1.5*u[i][1] - 0.5*u[i][2]);
				v[i][jend] = 2.0*v_bc - (1.5*v[i][1] - 0.5*v[i][2]);
				p[i][jend] = fluid_model->get_p_rhoT(rho[i][jend], T[i][jend]);
			}
			else{
				T[i][jend] = 1.5*T[i][njc] - 0.5*T[i][njc-1];
				rho[i][jend] = 1.5*rho[i][njc] - 0.5*rho[i][njc-1];
				u[i][jend] = 2.0*u_bc - (1.5*u[i][njc] - 0.5*u[i][njc-1]);
				v[i][jend] = 2.0*v_bc - (1.5*v[i][njc] - 0.5*v[i][njc-1]);
				p[i][jend] = fluid_model->get_p_rhoT(rho[i][jend], T[i][jend]);
			}
		}
	}

	if(face == global_left || face == global_right){
		if(face == global_left)
			iend = 0;
		else
			iend = nic + 1;

#pragma omp parallel for

		for(size_t j=start; j<=end; j++){
			if(face == global_left){
				T[iend][j] = 1.5*T[1][j] - 0.5*T[2][j];
				rho[iend][j] = 1.5*rho[1][j] - 0.5*rho[2][j];
				u[iend][j] = 2.0*u_bc - (1.5*u[1][j] - 0.5*u[2][j]);
				v[iend][j] = 2.0*v_bc - (1.5*v[1][j] - 0.5*v[2][j]);
				p[iend][j] = fluid_model->get_p_rhoT(rho[iend][j], T[iend][j]);
			}
			else{
				T[iend][j] = 1.5*T[nic][j] - 0.5*T[nic-1][j];
				rho[iend][j] = 1.5*rho[nic][j] - 0.5*rho[nic-1][j];
				u[iend][j] = 2.0*u_bc - (1.5*u[nic][j] - 0.5*u[nic-1][j]);
				v[iend][j] = 2.0*v_bc - (1.5*v[nic][j] - 0.5*v[nic-1][j]);
				p[iend][j] = fluid_model->get_p_rhoT(rho[iend][j], T[iend][j]);
			}
				
		}
	}
		
};



template<class Tx, class Tad>
BoundaryConditionWake<Tx, Tad>::BoundaryConditionWake(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
											 size_t val_face, size_t val_start, size_t val_end): BoundaryCondition<Tx, Tad>(){
	name = val_name;
	mesh = val_mesh;
	face = val_face;
	start = val_start;
	end = val_end;
		
	ni = mesh->ni;
	nj = mesh->nj;
	nic = ni - 1;
	njc = nj - 1;
};
	
template<class Tx, class Tad>
void BoundaryConditionWake<Tx, Tad>::apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
	size_t iend, jend;
	if(face == global_bottom || face == global_top){
		if(face == global_bottom)
			jend = 0;
		else
			jend = njc + 1;

#pragma omp parallel for
		for(size_t i=start; i<=end; i++){
			rho[i][jend] = rho[nic+1-i][1];
			u[i][jend] = u[nic+1-i][1];
			v[i][jend] = v[nic+1-i][1];
			p[i][jend] = p[nic+1-i][1];
			T[i][jend] = T[nic+1-i][1];
		}	

#pragma omp parallel for
		for(size_t i=start; i<=end; i++){
			rho[nic+1-i][jend] = rho[i][1];
			u[nic+1-i][jend] = u[i][1];
			v[nic+1-i][jend] = v[i][1];
			p[nic+1-i][jend] = p[i][1];
			T[nic+1-i][jend] = T[i][1];
		}
	}

	if(face == global_left || face == global_right){
		if(face == global_left)
			iend = 0;
		else
			iend = nic + 1;

#pragma omp parallel for
		for(size_t j=start; j<=end; j++){
			spdlog::get("console")->critical("Boundary condition not implemented!");
		}
	}
		
};

template<class Tx, class Tad>
BoundaryConditionOutflow<Tx, Tad>::BoundaryConditionOutflow(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> mesh, std::shared_ptr<Config<Tx>> config,std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
												   size_t val_face, size_t val_start, size_t val_end): BoundaryCondition<Tx, Tad>(){
	fluid_model = val_fluid_model;
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
	
template<class Tx, class Tad>
void BoundaryConditionOutflow<Tx, Tad>::apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
	size_t iend, jend;
	if(face == global_bottom || face == global_top){
		if(face == global_bottom)
			jend = 0;
		else
			jend = njc + 1;

		spdlog::get("console")->critical("Boundary not implemented!");
	}

	if(face == global_left || face == global_right){
		if(face == global_left)
			iend = 0;
		else
			iend = nic + 1;

#pragma omp parallel for
		for(size_t j=start; j<=end; j++){
			rho[iend][j] = rho[iend-1][j];
			u[iend][j] = u[iend-1][j];
			v[iend][j] = v[iend-1][j];
			p[iend][j] = p_inf;
			T[iend][j] = fluid_model->get_T_prho(p[iend][j], rho[iend][j]);
		}
	}
		
};



template<class Tx, class Tad>
BoundaryConditionPeriodic<Tx, Tad>::BoundaryConditionPeriodic(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> mesh, std::shared_ptr<Config<Tx>> config,std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
													 size_t val_face, size_t val_start, size_t val_end): BoundaryCondition<Tx, Tad>(){
	name = val_name;
	face = val_face;
	start = val_start;
	end = val_end;
	ni = mesh->ni;
	nj = mesh->nj;
	nic = ni - 1;
	njc = nj - 1;
};
	
template<class Tx, class Tad>
void BoundaryConditionPeriodic<Tx, Tad>::apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
	size_t iend, jend;
	if(face == global_bottom || face == global_top){
#pragma omp parallel for
		for(size_t i=start; i<=end; i++){
			rho[i][0] = rho[i][njc];
			u[i][0] = u[i][njc];
			v[i][0] = v[i][njc];
			p[i][0] = p[i][njc];
			T[i][0] = T[i][njc];
				
			rho[i][njc+1] = rho[i][1];
			u[i][njc+1] = u[i][1];
			v[i][njc+1] = v[i][1];
			p[i][njc+1] = p[i][1];
			T[i][njc+1] = T[i][1];
		}
	}

	if(face == global_left || face == global_right){
#pragma omp parallel for
		for(size_t j=start; j<=end; j++){
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
		}
	}
		
};

template<class Tx, class Tad>
BoundaryConditionIsothermalWall<Tx, Tad>::BoundaryConditionIsothermalWall(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> config,std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
																 size_t val_face, size_t val_start, size_t val_end, Tx val_u_bc, Tx val_v_bc, Tx val_T_bc): BoundaryCondition<Tx, Tad>(){
	fluid_model = val_fluid_model;
	name = val_name;
	mesh = val_mesh;
	face = val_face;
	start = val_start;
	end = val_end;
	ni = mesh->ni;
	nj = mesh->nj;
	nic = ni - 1;
	njc = nj - 1;

	u_bc = val_u_bc;
	v_bc = val_v_bc;
	T_bc = val_T_bc;
		
};
	
template<class Tx, class Tad>
void BoundaryConditionIsothermalWall<Tx, Tad>::apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
	size_t iend, jend;
	if(face == global_bottom || face == global_top){
		if(face == global_bottom)
			jend = 0;
		else
			jend = njc + 1;

#pragma omp parallel for
		for(size_t i=start; i<=end; i++){
			if(face == global_bottom){
				p[i][jend] = 1.5*p[i][1] - 0.5*p[i][2];
				u[i][jend] = 2.0*u_bc - (1.5*u[i][1] - 0.5*u[i][2]);
				v[i][jend] = 2.0*v_bc - (1.5*v[i][1] - 0.5*v[i][2]);
				T[i][jend] = T_bc;
				rho[i][jend] = fluid_model->get_rho_pT(p[i][jend], T[i][jend]);
			}
			else{
				p[i][jend] = 1.5*p[i][njc] - 0.5*p[i][njc-1];
				u[i][jend] = 2.0*u_bc - (1.5*u[i][njc] - 0.5*u[i][njc-1]);
				v[i][jend] = 2.0*v_bc - (1.5*v[i][njc] - 0.5*v[i][njc-1]);
				T[i][jend] = T_bc;
				rho[i][jend] = fluid_model->get_rho_pT(p[i][jend], T[i][jend]);
			}
		}
	}

	if(face == global_left || face == global_right){
		if(face == global_left)
			iend = 0;
		else
			iend = nic + 1;

#pragma omp parallel for
		for(size_t j=start; j<=end; j++){
			spdlog::get("console")->critical("Boundary condition not implemented!");
		}
	}
		
};

template<class Tx, class Tad>
void BoundaryContainer<Tx, Tad>::apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){
	for(auto&& bc : boundary_conditions)
		bc->apply(rho, u, v, p, T);
};

template<class Tx, class Tad>
size_t BoundaryContainer<Tx, Tad>::get_index(int idx, size_t face){
	if(idx >= 0){
		return idx;
	}
	else{
		if(face == global_bottom){
			return nic + 2 + idx;
		}
		else if(face == global_right){
			return njc + 2 + idx;
		}
		else if(face == global_top){
			return nic + 2 + idx;
		}
		else if(face == global_left){
			return njc + 2 + idx;
		}
		else{
			spdlog::get("console")->critical("something went wrong!");
		}
	}
};
	
template<class Tx, class Tad>
BoundaryContainer<Tx, Tad>::~BoundaryContainer(){
	for(auto&& bc : boundary_conditions)
		delete bc;
};
template<class Tx, class Tad>
BoundaryContainer<Tx, Tad>::BoundaryContainer(std::string filename, std::shared_ptr<Mesh<Tx, Tad>> mesh, std::shared_ptr<Config<Tx>> val_config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model){
	ni = mesh->ni;
	nj = mesh->nj;
	nic = ni - 1;
	njc = nj - 1;
	auto config = cpptoml::parse_file(filename);
	auto bcs = config->get_table_array("boundary");
	for (const auto& bc : *bcs){
		std::string name = bc->get_qualified_as<std::string>("name").value_or("boundary");
		std::string type = bc->get_qualified_as<std::string>("type").value_or("");
		std::string face = bc->get_qualified_as<std::string>("face").value_or("");
		size_t start = bc->get_qualified_as<int64_t>("start").value_or(0);
		size_t end = bc->get_qualified_as<int64_t>("end").value_or(0);

		size_t facei;
		if(face=="left") facei = global_left;
		if(face=="right") facei = global_right;
		if(face=="bottom") facei = global_bottom;
		if(face=="top") facei = global_top;
		end = get_index(end, facei);

		auto u_bc = bc->get_qualified_as<double>("u").value_or(0.0);
		auto v_bc = bc->get_qualified_as<double>("v").value_or(0.0);
		auto T_bc = bc->get_qualified_as<double>("T").value_or(0.0);
			
		//auto c = config->config;
		//spdlog::get("console")->debug("{}", c);
		//auto tmp_inf =  c->get_qualified_as<double>("freestream.rho_inf").value_or(1.0);
			
		if(type == "freestream"){
			auto boundarycondition = new BoundaryConditionFreestream<Tx, Tad>(name, mesh, val_config, val_fluid_model, facei, start, end);
			boundary_conditions.push_back(boundarycondition);
		}
		else if(type=="outflow"){
			auto boundarycondition = new BoundaryConditionOutflow<Tx, Tad>(name, mesh, val_config, val_fluid_model, facei, start, end);
			boundary_conditions.push_back(boundarycondition);
		}
		else if(type=="periodic"){
			auto boundarycondition = new BoundaryConditionPeriodic<Tx, Tad>(name, mesh, val_config, val_fluid_model, facei, start, end);
			boundary_conditions.push_back(boundarycondition);
		}
		else if(type == "slipwall"){
			auto boundarycondition = new BoundaryConditionInviscidWall<Tx, Tad>(name, mesh, val_config, val_fluid_model, facei, start, end);
			boundary_conditions.push_back(boundarycondition);
		}
		else if(type == "wall"){
			auto boundarycondition = new BoundaryConditionAdiabaticWall<Tx, Tad>(name, mesh, val_config, val_fluid_model, facei, start, end, u_bc, v_bc);
			boundary_conditions.push_back(boundarycondition);
		}
		else if(type == "isothermalwall"){
			auto boundarycondition = new BoundaryConditionIsothermalWall<Tx, Tad>(name, mesh, val_config, val_fluid_model, facei, start, end, u_bc, v_bc, T_bc);
			boundary_conditions.push_back(boundarycondition);
		}
		else if(type == "wake"){
			auto boundarycondition = new BoundaryConditionWake<Tx, Tad>(name, mesh, val_config, val_fluid_model, facei, start, end);
			boundary_conditions.push_back(boundarycondition);
		}
		else{
			spdlog::get("console")->info("Wrong type of BC.");
		}
	}
}

#if defined(ENABLE_ADOLC)
template class BoundaryContainer<double,adouble>;
template class BoundaryCondition<double,adouble>;
template class BoundaryConditionFreestream<double,adouble>;
template class BoundaryConditionOutflow<double,adouble>;
template class BoundaryConditionPeriodic<double,adouble>;
template class BoundaryConditionInviscidWall<double,adouble>;
template class BoundaryConditionAdiabaticWall<double,adouble>;
template class BoundaryConditionIsothermalWall<double,adouble>;
template class BoundaryConditionWake<double,adouble>;
#else
template class BoundaryContainer<double,double>;
template class BoundaryCondition<double,double>;
template class BoundaryConditionFreestream<double,double>;
template class BoundaryConditionOutflow<double,double>;
template class BoundaryConditionPeriodic<double,double>;
template class BoundaryConditionInviscidWall<double,double>;
template class BoundaryConditionAdiabaticWall<double,double>;
template class BoundaryConditionIsothermalWall<double,double>;
template class BoundaryConditionWake<double,double>;

template class BoundaryContainer<float,float>;
template class BoundaryCondition<float,float>;
template class BoundaryConditionFreestream<float,float>;
template class BoundaryConditionOutflow<float,float>;
template class BoundaryConditionPeriodic<float,float>;
template class BoundaryConditionInviscidWall<float,float>;
template class BoundaryConditionAdiabaticWall<float,float>;
template class BoundaryConditionIsothermalWall<float,float>;
template class BoundaryConditionWake<float,float>;
#endif
#endif
