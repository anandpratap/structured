template <class T, class Ti>
T*** allocate_3d_array(Ti nx, Ti ny, Ti nz){
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
void release_3d_array(T*** A, Ti nx, Ti ny, Ti nz){
    for (Ti i = 0; i < nx; ++i){
			for (Ti j = 0; j < ny; ++j){
				delete[] A[i][j];
			}
			delete[] A[i];
	}
    delete[] A;
}

template <class T, class Ti>
T** allocate_2d_array(Ti nx, Ti ny){
    T** A = new T*[nx];
    for(Ti i(0); i < nx; ++i){
		A[i] = new T[ny];
	}
    return A;
}

template <class T, class Ti>
void release_2d_array(T** A, Ti nx, Ti ny){
    for (Ti i = 0; i < nx; ++i){
		delete[] A[i];
	}
    delete[] A;
}


template <class T, class Ti>
T* allocate_1d_array(Ti nx){
    T *A = new T[nx];
    return A;
}

template <class T, class Ti>
	void release_1d_array(T* A, Ti nx){
    delete[] A;
}
