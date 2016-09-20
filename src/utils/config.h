#ifndef _CONFIG_H
#define _CONFIG_H
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
	Profiler(){
		timer_main = std::make_shared<Timer>();
		timer_residual = std::make_shared<Timer>();
		timer_jacobian = std::make_shared<Timer>();
		timer_linearsolver = std::make_shared<Timer>();
	}
	float current_time(){return timer_main->diff();}
	float update_time_residual(){t_residual += timer_residual->diff(); return timer_residual->diff();}
	void reset_time_residual(){timer_residual->reset();}

	float update_time_jacobian(){t_jacobian += timer_jacobian->diff();return timer_jacobian->diff();}
	void reset_time_jacobian(){timer_jacobian->reset();}

	float update_time_linearsolver(){t_linearsolver += timer_linearsolver->diff();return timer_linearsolver->diff();}
	void reset_time_linearsolver(){timer_linearsolver->reset();}
	void print(){
		float t_total = timer_main->diff();
		auto logger = spdlog::get("console");
		logger->info("----------");
		logger->info("Profiler");
		logger->info("total time = {}s", t_total);
		logger->info("time calculating residual = {}% ({}s)", t_residual/t_total*100, t_residual);
		logger->info("time calculating jacobian = {}% ({}s)", t_jacobian/t_total*100, t_jacobian);
		logger->info("time calculating linearsolver = {}% ({}s)", t_linearsolver/t_total*100, t_linearsolver);
	}
};

template<class Tx>
class ConfigFreestream{
public:
	Tx rho_inf, u_inf, v_inf, p_inf, mu_inf, pr_inf, T_inf, aoa;
	bool if_viscous;
	void set(std::shared_ptr<cpptoml::table> config){
		rho_inf =  config->get_qualified_as<double>("freestream.rho_inf").value_or(1.0);
		u_inf =  config->get_qualified_as<double>("freestream.u_inf").value_or(0.0);
		v_inf =  config->get_qualified_as<double>("freestream.v_inf").value_or(0.0);
		p_inf =  config->get_qualified_as<double>("freestream.p_inf").value_or(1.0/1.4);
		T_inf =  config->get_qualified_as<double>("freestream.T_inf").value_or(1.0/1.4);
		mu_inf =  config->get_qualified_as<double>("freestream.mu_inf").value_or(0.0);
		pr_inf =  config->get_qualified_as<double>("freestream.pr_inf").value_or(0.7);
		aoa =  config->get_qualified_as<double>("freestream.aoa").value_or(0.0)*M_PI/180.0;
		if_viscous = (mu_inf > 1e-15) ? true : false;
	};
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
	
	void set(std::shared_ptr<cpptoml::table> config){
		iteration_max = config->get_qualified_as<int64_t>("solver.iteration_max").value_or(1);

		order = config->get_qualified_as<int64_t>("solver.order").value_or(1);
		lhs_order = config->get_qualified_as<int64_t>("solver.lhs_order").value_or(order);
		
		time_accurate =  config->get_qualified_as<bool>("solver.time_accurate").value_or(false);
		scheme = config->get_qualified_as<std::string>("solver.scheme").value_or("forward_euler");

		flux = config->get_qualified_as<std::string>("solver.flux").value_or("ausm");

		tolerance = config->get_qualified_as<double>("solver.tolerance").value_or(1e10);

		cfl = config->get_qualified_as<double>("solver.cfl").value_or(1.0);
		cfl_ramp =  config->get_qualified_as<bool>("solver.cfl_ramp").value_or(false);
		cfl_ramp_iteration = config->get_qualified_as<int64_t>("solver.cfl_ramp_iteration").value_or(20);
		cfl_ramp_exponent = config->get_qualified_as<double>("solver.cfl_ramp_exponent").value_or(1.1);


		under_relaxation = config->get_qualified_as<double>("solver.under_relaxation").value_or(1.0);
		under_relaxation_ramp =  config->get_qualified_as<bool>("solver.under_relaxation_ramp").value_or(false);
		under_relaxation_ramp_iteration = config->get_qualified_as<int64_t>("solver.under_relaxation_ramp_iteration").value_or(20);
		under_relaxation_ramp_exponent = config->get_qualified_as<double>("solver.under_relaxation_ramp_exponent").value_or(1.1);

		dpdx = config->get_qualified_as<double>("source.dpdx").value_or(0.0);
		dpdy = config->get_qualified_as<double>("source.dpdy").value_or(0.0);
	};
};
template<class Tx>
class ConfigIO{
public:
	size_t stdout_frequency, fileout_frequency;
	bool restart;
	std::string label;
	void set(std::shared_ptr<cpptoml::table> config){
		restart =  config->get_qualified_as<bool>("io.restart").value_or(false);
		label = config->get_qualified_as<std::string>("io.label").value_or("flow");
		stdout_frequency = config->get_qualified_as<int64_t>("io.stdout_frequency").value_or(1);
		fileout_frequency = config->get_qualified_as<int64_t>("io.fileout_frequency").value_or(1);
	};
};

template<class Tx>
class ConfigGeometry{
public:
	size_t ni, nj, tail;
	std::string filename;
	std::string format;
	void set(std::shared_ptr<cpptoml::table> config){
		filename = config->get_qualified_as<std::string>("geometry.filename").value_or("grid.unf2");
		format = config->get_qualified_as<std::string>("geometry.format").value_or("grid.unf2");
		ni= config->get_qualified_as<int64_t>("geometry.ni").value_or(0);
		nj= config->get_qualified_as<int64_t>("geometry.nj").value_or(0);
		tail= config->get_qualified_as<int64_t>("geometry.tail").value_or(0);
	}
	~ConfigGeometry(){
	}
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
	
	Config(std::string val_filename, int val_argc, char *val_argv[]){
		argc = val_argc;
		argv = val_argv;
		freestream = std::make_shared<ConfigFreestream<Tx>>();
		solver = std::make_shared<ConfigSolver<Tx>>();
		io = std::make_shared<ConfigIO<Tx>>();
		geometry = std::make_shared<ConfigGeometry<Tx>>();

		filename = val_filename;
		config = cpptoml::parse_file(val_filename);
		freestream->set(config);
		solver->set(config);
		io->set(config);
		geometry->set(config);
		print();

		profiler = std::make_shared<Profiler>();
		
	};

	void print(){
		auto logger = spdlog::get("console");
		logger->info("Compiler: {}", __COMPILER__);
		logger->info("Compiled on {} at time {}", __DATE__, __TIME__);

		logger->info("---------------");
		
		logger->info("Configurations");
		PRINT_CONFIG(freestream->rho_inf);
		PRINT_CONFIG(freestream->u_inf);
		PRINT_CONFIG(freestream->v_inf);
		PRINT_CONFIG(freestream->p_inf);
		PRINT_CONFIG(freestream->T_inf);
		PRINT_CONFIG(freestream->mu_inf);
		PRINT_CONFIG(freestream->pr_inf);
		
		PRINT_CONFIG(solver->iteration_max);
		PRINT_CONFIG(solver->order);
		PRINT_CONFIG(solver->scheme);
		PRINT_CONFIG(solver->flux);
		PRINT_CONFIG(solver->time_accurate);
		PRINT_CONFIG(solver->tolerance);

		PRINT_CONFIG(solver->cfl);
		PRINT_CONFIG(solver->cfl_ramp);
		PRINT_CONFIG(solver->cfl_ramp_iteration);
		PRINT_CONFIG(solver->cfl_ramp_exponent);

		PRINT_CONFIG(solver->under_relaxation);
		PRINT_CONFIG(solver->under_relaxation_ramp);
		PRINT_CONFIG(solver->under_relaxation_ramp_iteration);
		PRINT_CONFIG(solver->under_relaxation_ramp_exponent);

		PRINT_CONFIG(io->restart);
		PRINT_CONFIG(io->label);
		PRINT_CONFIG(io->stdout_frequency);
		PRINT_CONFIG(io->fileout_frequency);

		PRINT_CONFIG(solver->dpdx);
		PRINT_CONFIG(solver->dpdy);
		
		PRINT_CONFIG(geometry->filename);
		PRINT_CONFIG(geometry->format);
		PRINT_CONFIG(geometry->ni);
		PRINT_CONFIG(geometry->nj);
		PRINT_CONFIG(geometry->tail);

		logger->info("---------------");

	};
};
#endif
