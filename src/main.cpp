#include "mesh.h"
#include "solver.h"
#include "solution.h"

int main(void){
	Mesh<double> m = Mesh<double>();
	std::cout<<m.ni<<std::endl;

	Mesh<double> mc = Mesh<double>(&m, 1, 1);
	std::cout<<mc.ni<<std::endl;


	Mesh<double> mc1 = Mesh<double>(&mc, 1, 1);
	std::cout<<mc1.ni<<std::endl;

	//Mesh<double> m = Mesh<double>();
	//std::cout<<m.ni<<std::endl;


	Solver<double> s = Solver<double>(&mc1);
	s.solve();
	write_solution(&mc1, "base.tec");
	
	return 0;
}

