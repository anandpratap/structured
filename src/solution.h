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
	T ***rhs;
	T **lhs;
	T ***q_aux;
	T *tmp_q;
	adouble *a_q, *a_rhs;
	arma::mat arma_rhs, arma_dq;
	arma::sp_mat arma_jac;
	
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
	rhs = allocate_3d_array<T>(nic, njc, nq);
	q_aux = allocate_3d_array<T>(nic, njc, naux);

	tmp_q = allocate_1d_array<T>(nic*njc*nq);
	a_q = allocate_1d_array<adouble>(nic*njc*nq);
	a_rhs = allocate_1d_array<adouble>(nic*njc*nq);
	
	arma_rhs = arma::mat(nic*njc*nq, 1);
	arma_dq = arma::mat(nic*njc*nq, 1);
	arma_jac = arma::sp_mat(nic*njc*nq, nic*njc*nq);

	
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			q[i][j][0] = 1.0;
			q[i][j][1] = 0.8;
			q[i][j][2] = 0.0;
			q[i][j][3] = 2.1057142857142863;
			for(uint k=0; k<nq; k++){
				rhs[i][j][k] = 0.0;
			}
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
	release_1d_array(tmp_q, nic*njc*nq);
	release_1d_array(a_q, nic*njc*nq);
	release_1d_array(a_rhs, nic*njc*nq);
	release_3d_array(rhs, nic, njc, nq);
	release_3d_array(q, nic, njc, nq);
	release_3d_array(q_aux, nic, njc, naux);
}


#endif
