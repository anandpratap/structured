#ifndef _DEF_CONFIG_H
#define _DEF_CONFIG_H
#include "common.h"
class Profiler{
public:
	std::shared_ptr<Timer> timer_main;
	std::shared_ptr<Timer> timer_residual;
	std::shared_ptr<Timer> timer_jacobian;
	std::shared_ptr<Timer> timer_linearsolver;
	float t_residual = 0.0f;
	float t_jacobian = 0.0f;
	float t_linearsolver = 0.0f;
	Profiler();
	float current_time();
	float update_time_residual();
	void reset_time_residual();

	float update_time_jacobian();
	void reset_time_jacobian();

	float update_time_linearsolver();
	void reset_time_linearsolver();
	void print();
};

template<class Tx>
class ConfigFreestream{
public:
	Tx rho_inf, u_inf, v_inf, p_inf, mu_inf, pr_inf, T_inf, aoa;
	bool if_viscous;
	void set(std::shared_ptr<cpptoml::table> config);
};
template<class Tx>
class ConfigSolver{
public:
	size_t order, cfl_ramp_iteration, under_relaxation_ramp_iteration, lhs_order;
	std::string scheme, flux;
	bool time_accurate, cfl_ramp, under_relaxation_ramp;
	Tx cfl, under_relaxation, cfl_ramp_exponent, under_relaxation_ramp_exponent;
	Tx tolerance;
	size_t iteration_max;
	Tx dpdx, dpdy;
	
	void set(std::shared_ptr<cpptoml::table> config);
};
template<class Tx>
class ConfigIO{
public:
	size_t stdout_frequency, fileout_frequency;
	bool restart;
	std::string label;
	void set(std::shared_ptr<cpptoml::table> config);
};

template<class Tx>
class ConfigGeometry{
public:
	size_t ni, nj, tail;
	std::string filename;
	std::string format;
	void set(std::shared_ptr<cpptoml::table> config);
};


template<class Tx>
class Config{
public:
	int argc;
	char **argv;
	std::shared_ptr<ConfigFreestream<Tx>> freestream;
	std::shared_ptr<ConfigSolver<Tx>> solver;
	std::shared_ptr<ConfigIO<Tx>> io;
	std::shared_ptr<ConfigGeometry<Tx>> geometry;
	std::shared_ptr<cpptoml::table> config;
	std::shared_ptr<Profiler> profiler;
	std::string filename;
	
	Config(std::string val_filename, int val_argc, char *val_argv[]);
	void print();
};
#endif
