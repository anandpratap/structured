#ifndef _UTILS_H
#define _UTILS_H

template <class T, class Ti>
T*** allocate_3d_array(const Ti nx, const Ti ny, const Ti nz){
    T*** A = new T**[nx];
    for(Ti i(0); i < nx; ++i){
			A[i] = new T*[ny];
			for(Ti j(0); j < ny; ++j){
				A[i][j] = new T[nz];					
				for(Ti k(0); k < nz; ++k){
					A[i][j][k]= 0.;
				}
			}
	}
    return A;
}
template <class T, class Ti>
void release_3d_array(T*** A, const Ti nx, const Ti ny, const Ti nz){
    for (Ti i = 0; i < nx; ++i){
			for (Ti j = 0; j < ny; ++j){
				delete[] A[i][j];
			}
			delete[] A[i];
	}
    delete[] A;
}

template <class T, class Ti>
T** allocate_2d_array(const Ti nx, const Ti ny){
    T** A = new T*[nx];
    for(Ti i(0); i < nx; ++i){
		A[i] = new T[ny];
	}
    return A;
}

template <class T, class Ti>
void release_2d_array(T** A, const Ti nx, const Ti ny){
    for (Ti i = 0; i < nx; ++i){
		delete[] A[i];
	}
    delete[] A;
}


template <class T, class Ti>
T* allocate_1d_array(const Ti nx){
    T *A = new T[nx];
    return A;
}

template <class T, class Ti>
	void release_1d_array(T* A, const Ti nx){
    delete[] A;
}


template<class Tad, class Ti>
void first_order_xi(const Ti ni, const Ti nj, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>& q, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>&  ql, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>&  qr){
	const auto njm = nj-1;
#pragma omp parallel for
	for(Ti i=0; i<ni; i++){
		for(Ti j=0; j<njm; j++){
			ql(i, j) = q(i, j+1);
			qr(i, j) = q(i+1, j+1);
		}
	}
}

template<class Tad, class Ti>
void first_order_eta(const Ti ni, const Ti nj, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>& q, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>&  ql, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>&  qr){
	const auto nim = ni-1;
#pragma omp parallel for
	for(Ti i=0; i<nim; i++){
		for(Ti j=0; j<nj; j++){
			ql(i, j) = q(i+1, j);
			qr(i, j) = q(i+1, j+1);
		}
	}
}


template<class Tad, class Ti>
void second_order_xi(const Ti ni, const Ti nj, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>& q, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>& ql, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>& qr){
	const auto njm = nj-1;

	ql.block(0, 0, ni, njm) = q.block(0, 1, ni, njm);
	qr.block(0, 0, ni, njm) = q.block(1, 1, ni, njm);
	const auto nim = ni-1;
	constexpr auto thm = 2.0/3.0;
	constexpr auto thp = 4.0/3.0;
	const auto eps = pow(10.0/nim, 3);
	static Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> f2  = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(ni, njm);
	static Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> a1 = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nim, njm);
	static Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> a2 = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nim, njm);
	static Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> f3qt = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nim, njm);

	f2.block(0, 0, ni, njm) = q.block(1, 1, ni, njm) - q.block(0, 1, ni, njm);
	a1.block(0, 0, nim, njm) = 3.0*f2.block(1, 0, nim, njm)*(f2.block(0,0,nim,njm));
	a2.block(0, 0, nim, njm) = 2.0*(f2.block(1, 0, nim, njm) - f2.block(0, 0, nim, njm))*(f2.block(1, 0, nim, njm) - f2.block(0, 0, nim, njm)) + a1.block(0, 0, nim, njm);

	f3qt = 0.25*(a1 + eps)/(a2 + eps);


	ql.block(1, 0, nim, njm) += f3qt*(thm*f2.block(0,0,nim,njm) + thp*f2.block(1,0,nim,njm));
	qr.block(0, 0, nim, njm) -= f3qt*(thp*f2.block(0,0,nim,njm) + thm*f2.block(1,0,nim,njm));
}

template<class Tad, class Ti>
void second_order_eta(const Ti ni, const Ti nj, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>& q, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>& ql, Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>& qr){
	const Ti nim = ni-1;

	ql.block(0, 0, nim, nj) = q.block(1, 0, nim, nj);
	qr.block(0, 0, nim, nj) = q.block(1, 1, nim, nj);
	const auto njm = nj-1;
	constexpr auto thm = 2.0/3.0;
	constexpr auto thp = 4.0/3.0;
	const auto eps = pow(10.0/njm, 3);
	static Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> f2  = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nim, nj);
	static Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> a1 = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nim, njm);
	static Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> a2 = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nim, njm);
	static Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic> f3qt = Eigen::Array<Tad, Eigen::Dynamic, Eigen::Dynamic>(nim, njm);


	f2.block(0, 0, nim, nj) = q.block(1, 1, nim, nj) - q.block(1, 0, nim, nj);

	a1.block(0, 0, nim, njm) = 3.0*f2.block(0, 1, nim, njm)*(f2.block(0,0,nim,njm));
	a2.block(0, 0, nim, njm) = 2.0*(f2.block(0, 1, nim, njm) - f2.block(0, 0, nim, njm))*(f2.block(0, 1, nim, njm) - f2.block(0, 0, nim, njm)) + a1.block(0, 0, nim, njm);
	f3qt = 0.25*(a1 + eps)/(a2 + eps);



	ql.block(0, 1, nim, njm) += f3qt*(thm*f2.block(0,0,nim,njm) + thp*f2.block(0,1,nim,njm));
	qr.block(0, 0, nim, njm) -= f3qt*(thp*f2.block(0,0,nim,njm) + thm*f2.block(0,1,nim,njm));
}


template<class T>
void primvars(const T Q[4], T *rho, T *u, T *v, T *p){
  const T tmp_rho = Q[0];
  const T tmp_u = Q[1]/tmp_rho;
  const T tmp_v = Q[2]/tmp_rho;
  *rho = Q[0];
  *u = Q[1]/Q[0];
  *v = Q[2]/Q[0];
  *p = (Q[3] - 0.5*tmp_rho*(tmp_u*tmp_u + tmp_v*tmp_v))*(GAMMA-1);
}


#endif
