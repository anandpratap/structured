#include "common.h"
#include "utils/mesh.h"
#include "solver/solver.h"
#include "solver/solution.h"
#include "utils/config.h"

int main(int argc, char *argv[]){
	spdlog::set_pattern("%v");

	auto logger = spdlog::stdout_logger_mt("console", true);
	logger->set_level(spdlog::level::debug);

	Config config("config.inp", argc, argv);
	
	Mesh<double> m = Mesh<double>(&config);
	Mesh<double> mc = Mesh<double>(&m, 1, 1);
	Mesh<double> mc1 = Mesh<double>(&mc, 1, 1);
	Mesh<double> mc2 = Mesh<double>(&mc1, 1, 1);

	Solver<double, adouble> s = Solver<double, adouble>(&mc, &config);
	s.solve();

	config.profiler->print();
	return 0;
}

