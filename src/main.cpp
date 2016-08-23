#include "common.h"
#include "utils/mesh.h"
#include "solver/solver.h"
#include "solver/solution.h"
#include "utils/config.h"

int main(int argc, char *argv[]){
	spdlog::set_pattern("%v");

	auto logger = spdlog::stdout_logger_mt("console", true);
	logger->set_level(spdlog::level::debug);

	auto config = std::make_shared<Config>("config.inp", argc, argv);
	auto m = std::make_shared<Mesh<double>>(config);
	auto mc = std::make_shared<Mesh<double>>(m, 1, 1);
	
#if defined(ENABLE_ADOLC)
	auto s = std::make_shared<Solver<double, adouble>>(m, config);
#else
	auto s = std::make_shared<Solver<double, double>>(m, config);
#endif
	s->solve();
	config->profiler->print();
	return 0;
}

