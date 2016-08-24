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


template<class T, class Ti>
void first_order_xi(const Ti ni, const Ti nj, Array2D<T>& q, Array2D<T>& ql, Array2D<T>& qr){
	const auto njm = nj-1;
#pragma omp parallel for
	for(Ti i=0; i<ni; i++){
		for(Ti j=0; j<njm; j++){
			ql(i, j) = q(i, j+1);
			qr(i, j) = q(i+1, j+1);
		}
	}
}

template<class T, class Ti>
void first_order_eta(const Ti ni, const Ti nj, Array2D<T>& q, Array2D<T>& ql, Array2D<T>& qr){
	const auto nim = ni-1;
#pragma omp parallel for
	for(Ti i=0; i<nim; i++){
		for(Ti j=0; j<nj; j++){
			ql(i, j) = q(i+1, j);
			qr(i, j) = q(i+1, j+1);
		}
	}
}


template<class T, class Ti>
void second_order_xi(const Ti ni, const Ti nj, Array2D<T>& q, Array2D<T>& ql, Array2D<T>& qr){
	const auto njm = nj-1;
#pragma omp parallel for
	for(Ti i=0; i<ni; i++){
		for(Ti j=0; j<njm; j++){
			ql(i, j) = q(i, j+1);
			qr(i, j) = q(i+1, j+1);
		}
	}
	const auto nim = ni-1;
	constexpr auto thm = 2.0/3.0;
	constexpr auto thp = 4.0/3.0;
	const auto eps = pow(10.0/nim, 3);
	static auto f2 = Array2D<T>(ni, njm);
	static auto a1 = Array2D<T>(nim, njm);
	static auto a2 = Array2D<T>(nim, njm);
	static auto f3qt = Array2D<T>(nim, njm);
#pragma omp parallel for
	for(Ti i=0; i<ni; i++){
		for(Ti j=0; j<njm; j++){
			f2(i, j) = q(i+1, j+1) - q(i, j+1);
		}
	}

#pragma omp parallel for
	for(Ti i=0; i<nim; i++){
		for(Ti j=0; j<njm; j++){
			a1(i, j) = 3.0*f2(i+1, j)*f2(i, j);
			a2(i, j) = 2*(f2(i+1, j) - f2(i, j))*(f2(i+1, j) - f2(i, j)) + a1(i, j);
			f3qt(i, j) = 0.25*(a1(i, j) + eps)/(a2(i, j) + eps);
		}
	}

#pragma omp parallel for
	for(Ti i=0; i<nim; i++){
		for(Ti j=0; j<njm; j++){
			ql(i+1, j) = ql(i+1, j) + f3qt(i, j)*(thm*f2(i, j) + thp*f2(i+1, j));
			qr(i, j) = qr(i, j) - f3qt(i, j)*(thp*f2(i, j) + thm*f2(i+1, j));
		}
	}
}

template<class T, class Ti>
void second_order_eta(const Ti ni, const Ti nj, Array2D<T>& q, Array2D<T>& ql, Array2D<T>& qr){
	const Ti nim = ni-1;

#pragma omp parallel for
	for(Ti i=0; i<nim; i++){
		for(Ti j=0; j<nj; j++){
			ql(i, j) = q(i+1, j);
			qr(i, j) = q(i+1, j+1);
		}
	}
	const auto njm = nj-1;
	constexpr auto thm = 2.0/3.0;
	constexpr auto thp = 4.0/3.0;
	const auto eps = pow(10.0/njm, 3);
	static auto f2 = Array2D<T>(nim, nj);
	static auto a1 = Array2D<T>(nim, njm);
	static auto a2 = Array2D<T>(nim, njm);
	static auto f3qt = Array2D<T>(nim, njm);

#pragma omp parallel for
	for(Ti i=0; i<nim; i++){
		for(Ti j=0; j<nj; j++){
			f2(i, j) = q(i+1, j+1) - q(i+1, j);
		}
	}

#pragma omp parallel for
	for(Ti i=0; i<nim; i++){
		for(Ti j=0; j<njm; j++){
			a1(i, j) = 3.0*f2(i, j+1)*f2(i, j);
			a2(i, j) = 2*(f2(i, j+1) - f2(i, j))*(f2(i, j+1) - f2(i, j)) + a1(i, j);
			f3qt(i, j) = 0.25*(a1(i, j) + eps)/(a2(i, j) + eps);
		}
	}
	
#pragma omp parallel for
	for(Ti i=0; i<nim; i++){
		for(Ti j=0; j<njm; j++){
			ql(i, j+1) = ql(i, j+1) + f3qt(i, j)*(thm*f2(i, j) + thp*f2(i, j+1));
			qr(i, j) = qr(i, j) - f3qt(i, j)*(thp*f2(i, j) + thm*f2(i, j+1));
		}
	}

}


template<class T>
void primvars(const T* Q, T *rho, T *u, T *v, T *p){
  const T tmp_rho = Q[0];
  const T tmp_u = Q[1]/tmp_rho;
  const T tmp_v = Q[2]/tmp_rho;
  *rho = Q[0];
  *u = Q[1]/Q[0];
  *v = Q[2]/Q[0];
  *p = (Q[3] - 0.5*tmp_rho*(tmp_u*tmp_u + tmp_v*tmp_v))*(GAMMA-1);
}


#endif
