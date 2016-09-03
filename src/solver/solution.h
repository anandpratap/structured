#ifndef _SOLUTION_H
#define _SOLUTION_H

#include "common.h"
#include "utils.h"
template<class T>
class Mesh;

#include "mesh.h"


template<class Tx>
class Solution{
 public:
	uint nic, njc, nq, naux, ntrans;
	Array3D<Tx> q;
	std::vector<std::string> q_name;
	Array3D<Tx> q_aux;
	std::vector<std::string> q_aux_name;
	Array2D<Tx> rho, u, v, p, T;
	std::vector<std::string> primvars_name;
	
 public:
	Solution(auto mesh);
	//	Solution(std::shared_ptr<Mesh<Tx>> mesh, std::shared_ptr<Mesh<Tx>> old_mesh, const uint nskipi=0, const uint nskipj=0, const uint refine=0);
	~Solution();
};

template<class Tx>
Solution<Tx>::Solution(auto mesh){
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

	for(uint i=0; i<naux; i++){
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
	
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			q[i][j][0] = 1.0;
			q[i][j][1] = 0.8;
			q[i][j][2] = 0.0;
			q[i][j][3] = 2.1057142857142863;
			for(uint tn=0; tn<ntrans; tn++){
				q[i][j][tn+4] = 1.0;
			}
		}
	}
}


//template<class Tx>
//Solution<Tx>::Solution(std::shared_ptr<Mesh<Tx>> mesh, std::shared_ptr<Mesh<Tx>> old_mesh, const uint nskipi, const uint nskipj, const uint refine): Solution<Tx>(mesh){
//	std::cout<<nq<<std::endl;
	// interpolate
//}


template<class Tx>
Solution<Tx>::~Solution(){
}


#endif
