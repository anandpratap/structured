#ifndef _SOLUTION_H
#define _SOLUTION_H

#include "common.h"
#include "utils.h"
template<class T>
class Mesh;

#include "mesh.h"


template<class T>
class Solution{
 public:
	uint nic, njc, nq, naux, ntrans;
	Array3D<T> q;
	Array3D<T> q_aux;
public:
	Solution(auto mesh);
	//	Solution(std::shared_ptr<Mesh<T>> mesh, std::shared_ptr<Mesh<T>> old_mesh, const uint nskipi=0, const uint nskipj=0, const uint refine=0);
	~Solution();
};

template<class T>
Solution<T>::Solution(auto mesh){
	nq = 4;
	ntrans = 0;
	naux = 10;
	nic = mesh->nic;
	njc = mesh->njc;

	q = Array3D<T>(nic, njc, nq+ntrans);
	q_aux = Array3D<T>(nic, njc, naux);

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


//template<class T>
//Solution<T>::Solution(std::shared_ptr<Mesh<T>> mesh, std::shared_ptr<Mesh<T>> old_mesh, const uint nskipi, const uint nskipj, const uint refine): Solution<T>(mesh){
//	std::cout<<nq<<std::endl;
	// interpolate
//}


template<class T>
Solution<T>::~Solution(){
}


#endif
