#ifndef _DEF_SOLVER_H
#define _DEF_SOLVER_H
#include "common.h"
#include "mesh.h"
#include "linearsolver.h"
#include "config.h"
#include "io.h"
#include "eulerequation.h"

template<class Tx>
void set_rarray(size_t size, Tx* __restrict__ dest, Tx* __restrict__ src);
template<class Tx>
void update_forward_euler(size_t size, Tx* __restrict__ q, Tx* __restrict__ rhs, Tx* __restrict__ dt);

template<class Tx, class To>
void update_rk4(size_t size, Tx* __restrict__ q_i, Tx* __restrict__ q, Tx* __restrict__ rhs, Tx* __restrict__ dt, To order);

template<class Tx, class Tad>
class Solver{
public:
	std::vector<std::shared_ptr<Mesh<Tx,Tad>>> mesh_list;
	Tx UNDER_RELAXATION;
	Tx CFL;
	std::shared_ptr<Config<Tx>> config;
	std::string label;
	std::shared_ptr<spdlog::logger> logger;
	std::shared_ptr<spdlog::logger> logger_convergence;	
	
	
	Solver(std::shared_ptr<Config<Tx>> config);
	~Solver();
	void add_mesh(std::shared_ptr<Mesh<Tx,Tad>> mesh);
	bool step(std::shared_ptr<Mesh<Tx,Tad>> mesh, size_t counter, Tx t);
	void solve();
};

#endif
