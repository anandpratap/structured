#include "common.h"
#include "utils/mesh.h"
#include "solver/solver.h"
#include "solver/solution.h"
#include "utils/config.h"

int main(int argc, char *argv[]){
#if defined(ENABLE_FLOAT)
	typedef float qtype;
#else
	typedef double qtype;
#endif

	spdlog::set_pattern("%v");

	auto logger = spdlog::stdout_logger_mt("console", true);
	logger->set_level(spdlog::level::debug);

	auto config = std::make_shared<Config<qtype>>("config.inp", argc, argv);
	auto m = std::make_shared<Mesh<qtype>>(config);

#if defined(ENABLE_ADOLC)
	auto s = std::make_shared<Solver<double, adouble>>(m, config);
#else
	auto s = std::make_shared<Solver<qtype, qtype>>(m, config);
#endif
	s->solve();
	config->profiler->print();
	return 0;
}

