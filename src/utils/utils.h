#ifndef _UTILS_H
#define _UTILS_H

template<class T>
void primvars(Array3D<T>& Q, Array2D<T>& rho, Array2D<T>& u, Array2D<T>& v, Array2D<T>& p, const uint shifti = 0, const uint shiftj = 0){
	uint nic = Q.extent(0);
	uint njc = Q.extent(1);
#pragma omp parallel for
	for(uint i=0; i<nic; i++){
		for(uint j=0; j<njc; j++){
			T tmp_rho, tmp_u, tmp_v;
			tmp_rho = Q[i][j][0];
			tmp_u = Q[i][j][1]/tmp_rho;
			tmp_v = Q[i][j][2]/tmp_rho;
			rho[i+shifti][j+shiftj] = tmp_rho;
			u[i+shifti][j+shiftj] = tmp_u;
			v[i+shifti][j+shiftj] = tmp_v;
			p[i+shifti][j+shiftj] = (Q[i][j][3] - 0.5*tmp_rho*(tmp_u*tmp_u + tmp_v*tmp_v))*(GAMMA-1);
		}
	}
  
}


#endif
