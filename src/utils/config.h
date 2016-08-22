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

	void update_time_residual(){t_residual += timer_residual->diff();}
	void reset_time_residual(){timer_residual->reset();}

	void update_time_jacobian(){t_jacobian += timer_jacobian->diff();}
	void reset_time_jacobian(){timer_jacobian->reset();}

	void update_time_linearsolver(){t_linearsolver += timer_linearsolver->diff();}
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

class ConfigFreestream{
 public:
	double rho_inf, u_inf, v_inf, p_inf;
	void set(std::shared_ptr<cpptoml::table> config){
		rho_inf =  config->get_qualified_as<double>("freestream.rho_inf").value_or(1.0);
		u_inf =  config->get_qualified_as<double>("freestream.u_inf").value_or(0.0);
		v_inf =  config->get_qualified_as<double>("freestream.v_inf").value_or(0.0);
		p_inf =  config->get_qualified_as<double>("freestream.p_inf").value_or(1.0/1.4);
	};
};

class ConfigSolver{
 public:
	int order, cfl_ramp_iteration, under_relaxation_ramp_iteration;
	std::string scheme;
	bool time_accurate, cfl_ramp, under_relaxation_ramp;
	double cfl, under_relaxation, cfl_ramp_exponent, under_relaxation_ramp_exponent;
	double tolerance;
	
	void set(std::shared_ptr<cpptoml::table> config){
		order = config->get_qualified_as<int64_t>("solver.order").value_or(1);
		time_accurate =  config->get_qualified_as<bool>("solver.time_accurate").value_or(false);
		scheme = config->get_qualified_as<std::string>("solver.scheme").value_or("forward_euler");

		tolerance = config->get_qualified_as<double>("solver.tolerance").value_or(1e10);

		cfl = config->get_qualified_as<double>("solver.cfl").value_or(1.0);
		cfl_ramp =  config->get_qualified_as<bool>("solver.cfl_ramp").value_or(false);
		cfl_ramp_iteration = config->get_qualified_as<int64_t>("solver.cfl_ramp_iteration").value_or(20);
		cfl_ramp_exponent = config->get_qualified_as<double>("solver.cfl_ramp_exponent").value_or(1.1);


		under_relaxation = config->get_qualified_as<double>("solver.under_relaxation").value_or(1.0);
		under_relaxation_ramp =  config->get_qualified_as<bool>("solver.under_relaxation_ramp").value_or(false);
		under_relaxation_ramp_iteration = config->get_qualified_as<int64_t>("solver.under_relaxation_ramp_iteration").value_or(20);
		under_relaxation_ramp_exponent = config->get_qualified_as<double>("solver.under_relaxation_ramp_exponent").value_or(1.1);
		
	};
};

class ConfigIO{
 public:
	int stdout_frequency, fileout_frequency;
	bool restart;
	std::string label;
	void set(std::shared_ptr<cpptoml::table> config){
		restart =  config->get_qualified_as<bool>("io.restart").value_or(false);
		label = config->get_qualified_as<std::string>("io.label").value_or("flow");
		stdout_frequency = config->get_qualified_as<int64_t>("io.stdout_frequency").value_or(1);
		fileout_frequency = config->get_qualified_as<int64_t>("io.fileout_frequency").value_or(1);
	};
};


class ConfigGeometry{
 public:
	int ni, nj, tail;
	std::string filename;
	void set(std::shared_ptr<cpptoml::table> config){
		filename = config->get_qualified_as<std::string>("geometry.filename").value_or("grid.unf2");
		ni= config->get_qualified_as<int64_t>("geometry.ni").value_or(0);
		nj= config->get_qualified_as<int64_t>("geometry.nj").value_or(0);
		tail= config->get_qualified_as<int64_t>("geometry.tail").value_or(0);
	}
};

class Config{
 public:
	int argc;
	char **argv;
	std::shared_ptr<ConfigFreestream> freestream;
	std::shared_ptr<ConfigSolver> solver;
	std::shared_ptr<ConfigIO> io;
	std::shared_ptr<ConfigGeometry> geometry;
	std::shared_ptr<cpptoml::table> config;
	std::shared_ptr<Profiler> profiler;
	
	Config(std::string filename, int val_argc, char *val_argv[]){
		argc = val_argc;
		argv = val_argv;
		freestream = std::make_shared<ConfigFreestream>();
		solver = std::make_shared<ConfigSolver>();
		io = std::make_shared<ConfigIO>();
		geometry = std::make_shared<ConfigGeometry>();
		
		config = cpptoml::parse_file(filename);
		freestream->set(config);
		solver->set(config);
		io->set(config);
		geometry->set(config);
		print();

		profiler = std::make_shared<Profiler>();
		
	};

	void print(){
		auto logger = spdlog::get("console");

		logger->info("---------------");
		logger->info("Freestream configuration");
		logger->info("(rho, u, v, p) = ({}, {}, {}, {})", freestream->rho_inf, freestream->u_inf,
					 freestream->v_inf, freestream->p_inf);
		logger->info("---------------");

		logger->info("---------------");
		logger->info("Solver configuration");
		logger->info("order: {}", solver->order);
		logger->info("scheme: {}", solver->scheme);
		logger->info("time_accurate: {}", solver->time_accurate);
		logger->info("tolerance: {}", solver->tolerance);

		logger->info("cfl: {}", solver->cfl);
		logger->info("cfl_ramp: {}", solver->cfl_ramp);
		logger->info("cfl_ramp_iteration: {}", solver->cfl_ramp_iteration);
		logger->info("cfl_ramp_exponent: {}", solver->cfl_ramp_exponent);

		logger->info("under_relaxation: {}", solver->under_relaxation);
		logger->info("under_relaxation_ramp: {}", solver->under_relaxation_ramp);
		logger->info("under_relaxation_ramp_iteration: {}", solver->under_relaxation_ramp_iteration);
		logger->info("under_relaxation_ramp_exponent: {}", solver->under_relaxation_ramp_exponent);

		logger->info("---------------");


		logger->info("---------------");
		logger->info("IO configuration");
		logger->info("label = {}", io->label);
		logger->info("---------------");

		logger->info("---------------");
		logger->info("Geometry configuration");
		logger->info("filename = {}", geometry->filename);
		logger->info("ni = {}", geometry->ni);
		logger->info("nj = {}", geometry->nj);
		logger->info("tail = {}", geometry->tail);

		logger->info("---------------");

	};
};


#endif
