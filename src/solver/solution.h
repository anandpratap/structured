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
	uint nic, njc, nq, naux;
	T ***q;
	T ***q_aux;
public:
	Solution(const Mesh<T>* mesh);
	Solution(const Mesh<T>* mesh, const Mesh<T>* old_mesh, const uint nskipi=0, const uint nskipj=0, const uint refine=0);
	~Solution();
};

template<class T>
Solution<T>::Solution(const Mesh<T>* mesh){
	nq = 4;
	naux = 10;
	nic = mesh->nic;
	njc = mesh->njc;

	q = allocate_3d_array<T>(nic, njc, nq);
	q_aux = allocate_3d_array<T>(nic, njc, naux);

	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			q[i][j][0] = 1.0;
			q[i][j][1] = 0.8;
			q[i][j][2] = 0.0;
			q[i][j][3] = 2.1057142857142863;
		}
	}
}


template<class T>
Solution<T>::Solution(const Mesh<T>* mesh, const Mesh<T>* old_mesh, const uint nskipi, const uint nskipj, const uint refine): Solution<T>(mesh){
	std::cout<<nq<<std::endl;
	// interpolate
}


template<class T>
Solution<T>::~Solution(){
	release_3d_array(q, nic, njc, nq);
	release_3d_array(q_aux, nic, njc, naux);
}


#endif
