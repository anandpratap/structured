#include "common.h"
#include "def_mesh.h"
#include "def_solver.h"
#include "def_solution.h"
#include "def_config.h"
#include <fenv.h>

int main(int argc, char *argv[]){
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	
#if defined(ENABLE_FLOAT)
	using qtype = float;
#else
	using qtype = double;
#endif

#if defined(ENABLE_ADOLC)
	using aqtype = adouble;
#else
	using aqtype = qtype;
#endif

	
	cmdline::parser parser;
	spdlog::set_pattern("%v");
	parser.add<std::string>("config", 'c', "configuration file name", false, "config.inp");
	parser.add("help", '\0', "help");
	parser.parse(argc, argv);
	if(parser.exist("help")){std::cout << parser.usage(); exit(0);}

	auto logger = spdlog::stdout_logger_mt("console", true);
	logger->set_level(spdlog::level::debug);

	auto config = std::make_shared<Config<qtype>>(parser.get<std::string>("config"), argc, argv);
	auto m = std::make_shared<Mesh<qtype, aqtype>>(config);
	m->label = "";
	m->setup();
	auto s = std::make_shared<Solver<qtype, aqtype>>(config);

	s->add_mesh(m);
	s->solve();
	config->profiler->print();
	return 0;
}

