#ifndef _ADJOINT_H
#define _ADJOINT_H
#include "common.h"
#include "mesh.h"

template<class Tx, class Tad>
class AdjointSolver{
public:
	AdjointSolver();
	~AdjointSolver();
	void solve(std::shared_ptr<Mesh<Tx,Tad>> mesh);
	void calc_sensitivity(std::shared_ptr<Mesh<Tx,Tad>> mesh);
};

#endif
