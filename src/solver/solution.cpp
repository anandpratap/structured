#ifndef _SOLUTION_CPP
#define _SOLUTION_CPP
#include "solution.h"
#include "common.h"
#include "mesh.h"
template<class Tx, class Tad>
Solution<Tx, Tad>::Solution(std::shared_ptr<Mesh<Tx,Tad>> mesh){
	nq = 4;
	ntrans = 0;
	naux = 10;

	q_name = std::vector<std::string>(nq);
	q_aux_name = std::vector<std::string>(naux);
	primvars_name = std::vector<std::string>(5);
	
	q_name[0] = "rho";
	q_name[1] = "rhou";
	q_name[2] = "rhov";
	q_name[3] = "rhoE";

	for(size_t i=0; i<naux; i++){
		q_aux_name[i] = "q_aux_" + std::to_string(i);
	}
	
	nic = mesh->nic;
	njc = mesh->njc;

	q = Array3D<Tx>(nic, njc, nq+ntrans);
	q_aux = Array3D<Tx>(nic, njc, naux);

	rho = Array2D<Tx>(nic, njc);
	u = Array2D<Tx>(nic, njc);
	v = Array2D<Tx>(nic, njc);
	p = Array2D<Tx>(nic, njc);
	T = Array2D<Tx>(nic, njc);
	
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			q[i][j][0] = 1.0;
			q[i][j][1] = 0.8;
			q[i][j][2] = 0.0;
			q[i][j][3] = 2.1057142857142863;
			for(size_t tn=0; tn<ntrans; tn++){
				q[i][j][tn+4] = 1.0;
			}
		}
	}


	nt = nic*njc*(nq+ntrans);
	dt = Array3D<Tx>(nic, njc, nq+ntrans);
	rhs = Array3D<Tx>(nic, njc, nq+ntrans);
	q_tmp = Array3D<Tx>(nic, njc, nq+ntrans);
	a_q = Array3D<Tad>(nic, njc, nq+ntrans);
	a_rhs = Array3D<Tad>(nic, njc, nq+ntrans);
}


//template<class Tx, class Tad>
//Solution<Tx>::Solution(std::shared_ptr<Mesh<Tx>> mesh, std::shared_ptr<Mesh<Tx>> old_mesh, const size_t nskipi, const size_t nskipj, const size_t refine): Solution<Tx>(mesh){
//	std::cout<<nq<<std::endl;
// interpolate
//}


template<class Tx, class Tad>
Solution<Tx, Tad>::~Solution(){
}

#if defined(ENABLE_ADOLC)
template class Solution<double, adouble>;
#else
template class Solution<double, double>;
template class Solution<float, float>;
#endif
#endif
