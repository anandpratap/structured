#ifndef _SOLUTION_H
#define _SOLUTION_H
#include "common.h"
#include "mesh.h"
template<class Tx, class Tad>
class Solution{
public:
	size_t nic, njc, nq, naux, ntrans;
	Array3D<Tx> q;
	std::vector<std::string> q_name;
	Array3D<Tx> q_aux;
	std::vector<std::string> q_aux_name;
	Array2D<Tx> rho, u, v, p, T;
	std::vector<std::string> primvars_name;

	int nnz;
	int repeat = 0;
	unsigned int *rind = nullptr;
	unsigned int *cind = nullptr;
	double *values = nullptr;
	int options[4] = {0,0,0,0};

	size_t nt;
	Array3D<Tx> rhs;
	Tx **lhs;
	Array3D<Tx> dt;
	Array3D<Tx> q_tmp;
	Array3D<Tad> a_q, a_rhs;
public:
	Solution(std::shared_ptr<Mesh<Tx,Tad>> mesh);
	//	Solution(std::shared_ptr<Mesh<Tx>> mesh, std::shared_ptr<Mesh<Tx>> old_mesh, const size_t nskipi=0, const size_t nskipj=0, const size_t refine=0);
	~Solution();
};
#endif
