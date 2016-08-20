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

	Mesh<double> mc2 = Mesh<double>(&mc1, 1, 1);
	std::cout<<mc2.ni<<std::endl;


	//Mesh<double> m = Mesh<double>();
	//std::cout<<m.ni<<std::endl;


	Solver<double> s = Solver<double>(&mc);
	s.solve();
	write_solution(&mc, "base.tec");
	
	return 0;
}

