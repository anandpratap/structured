#ifndef _CONFIG_CPP
#define _CONFIG_CPP
#include "common.h"
#include "config.h"
Profiler::Profiler(){
		timer_main = std::make_shared<Timer>();
		timer_residual = std::make_shared<Timer>();
		timer_jacobian = std::make_shared<Timer>();
		timer_linearsolver = std::make_shared<Timer>();
}
float Profiler::current_time(){return timer_main->diff();}
float Profiler::update_time_residual(){t_residual += timer_residual->diff(); return timer_residual->diff();}
void Profiler::reset_time_residual(){timer_residual->reset();}

float Profiler::update_time_jacobian(){t_jacobian += timer_jacobian->diff();return timer_jacobian->diff();}
void Profiler::reset_time_jacobian(){timer_jacobian->reset();}

float Profiler::update_time_linearsolver(){t_linearsolver += timer_linearsolver->diff();return timer_linearsolver->diff();}
void Profiler::reset_time_linearsolver(){timer_linearsolver->reset();}
void Profiler::print(){
		float t_total = timer_main->diff();
		auto logger = spdlog::get("console");
		logger->info("----------");
		logger->info("Profiler");
		logger->info("total time = {}s", t_total);
		logger->info("time calculating residual = {}% ({}s)", t_residual/t_total*100, t_residual);
		logger->info("time calculating jacobian = {}% ({}s)", t_jacobian/t_total*100, t_jacobian);
		logger->info("time calculating linearsolver = {}% ({}s)", t_linearsolver/t_total*100, t_linearsolver);
}

template<class Tx>
void ConfigFreestream<Tx>::set(std::shared_ptr<cpptoml::table> config){
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

template<class Tx>
void ConfigDesign<Tx>::set(std::shared_ptr<cpptoml::table> config){
	perturb_mode = config->get_qualified_as<std::string>("design.perturb_mode").value_or("none");
}

template<class Tx>
void ConfigSolver<Tx>::set(std::shared_ptr<cpptoml::table> config){
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
template<class Tx>
void ConfigIO<Tx>::set(std::shared_ptr<cpptoml::table> config){
		restart =  config->get_qualified_as<bool>("io.restart").value_or(false);
		label = config->get_qualified_as<std::string>("io.label").value_or("flow");
		stdout_frequency = config->get_qualified_as<int64_t>("io.stdout_frequency").value_or(1);
		fileout_frequency = config->get_qualified_as<int64_t>("io.fileout_frequency").value_or(1);
};

template<class Tx>
void ConfigGeometry<Tx>::set(std::shared_ptr<cpptoml::table> config){
		filename = config->get_qualified_as<std::string>("geometry.filename").value_or("grid.unf2");
		format = config->get_qualified_as<std::string>("geometry.format").value_or("grid.unf2");
		ni= config->get_qualified_as<int64_t>("geometry.ni").value_or(0);
		nj= config->get_qualified_as<int64_t>("geometry.nj").value_or(0);
		tail= config->get_qualified_as<int64_t>("geometry.tail").value_or(0);
}


template<class Tx>
Config<Tx>::Config(std::string val_filename, int val_argc, char *val_argv[]){
		argc = val_argc;
		argv = val_argv;
		freestream = std::make_shared<ConfigFreestream<Tx>>();
		design = std::make_shared<ConfigDesign<Tx>>();
		solver = std::make_shared<ConfigSolver<Tx>>();
		io = std::make_shared<ConfigIO<Tx>>();
		geometry = std::make_shared<ConfigGeometry<Tx>>();

		filename = val_filename;
		config = cpptoml::parse_file(val_filename);
		freestream->set(config);
		design->set(config);
		solver->set(config);
		io->set(config);
		geometry->set(config);
		print();

		profiler = std::make_shared<Profiler>();
		
};
template<class Tx>
void Config<Tx>::print(){
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

		PRINT_CONFIG(design->perturb_mode);
		logger->info("---------------");

};

template class Config<double>;
template class ConfigFreestream<double>;
template class ConfigIO<double>;
template class ConfigGeometry<double>;
template class ConfigSolver<double>;


template class Config<float>;
template class ConfigFreestream<float>;
template class ConfigIO<float>;
template class ConfigGeometry<float>;
template class ConfigSolver<float>;
#endif
