#ifndef _DEF_BC_H
#define _DEF_BC_H
#include "common.h"
#include "def_config.h"
#include "def_fluid.h"
#include "def_mesh.h"

#define global_bottom 0
#define global_right 1
#define global_top 2
#define global_left 3

template<class Tx, class Tad>
class BoundaryCondition{
public:
	virtual void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T){};
};

template<class Tx, class Tad>
class BoundaryConditionFreestream: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	size_t face;
	size_t start, end;
	size_t ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	std::shared_ptr<FluidModel<Tx, Tad>> fluid_model;
	BoundaryConditionFreestream(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> mesh, std::shared_ptr<Config<Tx>> config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
								size_t val_face, size_t val_start, size_t val_end);
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T);
};

template<class Tx, class Tad>
class BoundaryConditionInviscidWall: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	size_t face;
	size_t start, end;
	size_t ni, nj, nic, njc;
	std::shared_ptr<Mesh<Tx, Tad>> mesh;
	std::shared_ptr<FluidModel<Tx, Tad>> fluid_model;
	BoundaryConditionInviscidWall(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
								  size_t val_face, size_t val_start, size_t val_end);	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T);
};

template<class Tx, class Tad>
class BoundaryConditionAdiabaticWall: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	size_t face;
	size_t start, end;
	size_t ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	Tx u_bc, v_bc;
	std::shared_ptr<Mesh<Tx, Tad>> mesh;
	std::shared_ptr<FluidModel<Tx, Tad>> fluid_model;
	BoundaryConditionAdiabaticWall(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
								   size_t val_face, size_t val_start, size_t val_end, Tx val_u_bc, Tx val_v_bc);	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T);
};


template<class Tx, class Tad>
class BoundaryConditionWake: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	size_t face;
	size_t start, end;
	size_t ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	std::shared_ptr<Mesh<Tx, Tad>> mesh;
	BoundaryConditionWake(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
						  size_t val_face, size_t val_start, size_t val_end);	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T);
};

template<class Tx, class Tad>
class BoundaryConditionOutflow: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	size_t face;
	size_t start, end;
	size_t ni, nj, nic, njc;
	Tx rho_inf, u_inf, v_inf, p_inf;
	std::shared_ptr<FluidModel<Tx, Tad>> fluid_model;
	BoundaryConditionOutflow(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> mesh, std::shared_ptr<Config<Tx>> config,std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
							 size_t val_face, size_t val_start, size_t val_end);	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T);
};



template<class Tx, class Tad>
class BoundaryConditionPeriodic: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	size_t face;
	size_t start, end;
	size_t ni, nj, nic, njc;
	BoundaryConditionPeriodic(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> mesh, std::shared_ptr<Config<Tx>> config,std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
							  size_t val_face, size_t val_start, size_t val_end);	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T);
};


template<class Tx, class Tad>
class BoundaryConditionIsothermalWall: public BoundaryCondition<Tx, Tad>{
public:
	std::string name;
	size_t face;
	size_t start, end;
	size_t ni, nj, nic, njc;
	Tx u_bc, v_bc, T_bc;
	std::shared_ptr<Mesh<Tx, Tad>> mesh;
	std::shared_ptr<FluidModel<Tx, Tad>> fluid_model;
	BoundaryConditionIsothermalWall(std::string val_name, std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> config,std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model,
									size_t val_face, size_t val_start, size_t val_end, Tx val_u_bc, Tx val_v_bc, Tx val_T_bc);	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T);
};

template<class Tx, class Tad>
class BoundaryContainer{
public:
	size_t nic, njc, ni, nj;
	std::vector<BoundaryCondition<Tx, Tad>*> boundary_conditions;
	BoundaryContainer(std::string filename, std::shared_ptr<Mesh<Tx, Tad>> mesh, std::shared_ptr<Config<Tx>> val_config, std::shared_ptr<FluidModel<Tx, Tad>> val_fluid_model);	
	~BoundaryContainer();
	
	void apply(Array2D<Tad>& rho, Array2D<Tad>& u, Array2D<Tad>& v, Array2D<Tad>& p, Array2D<Tad>& T);
	size_t get_index(int idx, size_t face);	
};
#endif
