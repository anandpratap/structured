#ifndef _IO_H
#define _IO_H
#include "common.h"
#include "mesh.h"
#include "config.h"
#include "solution.h"
#include "eulerequation.h"

template<class Tq>
Tq value(const Tq x);

#if defined(ENABLE_ADOLC)
double value(const adouble x);
#endif

template<class Tx, class Tad>
class IOManager{
public:
	std::shared_ptr<Config<Tx>> config;
	std::shared_ptr<Mesh<Tx, Tad>> mesh;
	std::string label;
	IOManager(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config);
	~IOManager();
	void write(const size_t iteration);	
	void write_tecplot();	
	void write_tecplot_adjoint();	
	void write_npz();
	void write_restart();
	void read_restart();
	void write_surface();
};

#endif
