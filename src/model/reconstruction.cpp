#ifndef __RECONSTRUCTION__CPP
#define __RECONSTRUCTION__CPP
#include "common.h"
#include "reconstruction.h"
template<class Tx, class Tad>
void Reconstruction<Tx, Tad>::evaluate_chi(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){};
template<class Tx, class Tad>
void Reconstruction<Tx, Tad>::evaluate_eta(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){};
/*!
  \brief First order reconstruction

  Simple first order reconstruction where left and right states are taken to be the values of 
  left and right cell.
*/
template<class Tx, class Tad>
ReconstructionFirstOrder<Tx, Tad>::ReconstructionFirstOrder(const size_t val_ni, const size_t val_nj): Reconstruction<Tx, Tad>(){
	ni = val_ni;
	nj = val_nj;
	nic = ni - 1;
	njc = nj - 1;
};
/*!@copydoc
  \f{eqnarray*}{
  q_{left, i, j} = q_{i, j+1} \\
  q_{right, i, j} = q_{i+1, j+1} \\
  \f} 
*/	
template<class Tx, class Tad>
void ReconstructionFirstOrder<Tx, Tad>::evaluate_chi(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
#pragma omp parallel for
	for(size_t i=0; i<ni; i++){
		for(size_t j=0; j<njc; j++){
			ql[i][j] = q[i][j+1];
			qr[i][j] = q[i+1][j+1];
		}
	}
};
	
template<class Tx, class Tad>
void ReconstructionFirstOrder<Tx, Tad>::evaluate_eta(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
#pragma omp parallel for
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<nj; j++){
			ql[i][j] = q[i+1][j];
			qr[i][j] = q[i+1][j+1];
		}
	}
};

/*!
  \brief Second order reconstruction

  Simple second order reconstruction.
*/
template<class Tx, class Tad>
ReconstructionSecondOrder<Tx, Tad>::ReconstructionSecondOrder(const size_t val_ni, const size_t val_nj): Reconstruction<Tx, Tad>(){
	ni = val_ni;
	nj = val_nj;
	nic = ni - 1;
	njc = nj - 1;

	eps_chi = pow(10.0/nic, 3);
	eps_eta = pow(10.0/njc, 3); 
		
	f2_chi = Array2D<Tad>(ni, njc);
	a1_chi = Array2D<Tad>(nic, njc);
	a2_chi = Array2D<Tad>(nic, njc);
	f3qt_chi = Array2D<Tad>(nic, njc);


	f2_eta = Array2D<Tad>(nic, nj);
	a1_eta = Array2D<Tad>(nic, njc);
	a2_eta = Array2D<Tad>(nic, njc);
	f3qt_eta = Array2D<Tad>(nic, njc);
		
};
template<class Tx, class Tad>
void ReconstructionSecondOrder<Tx, Tad>::evaluate_chi(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
	auto f2 = f2_chi;
	auto a1 = a1_chi;
	auto a2 = a2_chi;
	auto f3qt = f3qt_chi;
	auto eps = eps_chi;
#pragma omp parallel  for
	for(size_t i=0; i<ni; i++){
		for(size_t j=0; j<njc; j++){
			ql[i][j] = q[i][j+1];
			qr[i][j] = q[i+1][j+1];
		}
	}
#pragma omp parallel  for
	for(size_t i=0; i<ni; i++){
		for(size_t j=0; j<njc; j++){
			f2[i][j] = q[i+1][j+1] - q[i][j+1];
		}
	}
		
#pragma omp parallel  for
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			a1[i][j] = 3.0*f2[i+1][j]*f2[i][j];
			a2[i][j] = 2.0*(f2[i+1][j] - f2[i][j])*(f2[i+1][j] - f2[i][j]) + a1[i][j];
			f3qt[i][j] = 0.25*(a1[i][j] + eps)/(a2[i][j] + eps);
		}
	}
		
#pragma omp parallel  for
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			ql[i+1][j] = ql[i+1][j] + f3qt[i][j]*(thm*f2[i][j] + thp*f2[i+1][j]);
			qr[i][j] = qr[i][j] - f3qt[i][j]*(thp*f2[i][j] + thm*f2[i+1][j]);
		}
	}
};

template<class Tx, class Tad>
void ReconstructionSecondOrder<Tx, Tad>::evaluate_eta(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
	auto f2 = f2_eta;
	auto a1 = a1_eta;
	auto a2 = a2_eta;
	auto f3qt = f3qt_eta;
	auto eps = eps_eta;
#pragma omp parallel  for
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<nj; j++){
			ql[i][j] = q[i+1][j];
			qr[i][j] = q[i+1][j+1];
		}
	}
#pragma omp parallel  for
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<nj; j++){
			f2[i][j] = q[i+1][j+1] - q[i+1][j];
		}
	}
		
#pragma omp parallel  for
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			a1[i][j] = 3.0*f2[i][j+1]*f2[i][j];
			a2[i][j] = 2.0*(f2[i][j+1] - f2[i][j])*(f2[i][j+1] - f2[i][j]) + a1[i][j];
			f3qt[i][j] = 0.25*(a1[i][j] + eps)/(a2[i][j] + eps);
		}
	}
		
#pragma omp parallel  for
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			ql[i][j+1] = ql[i][j+1] + f3qt[i][j]*(thm*f2[i][j] + thp*f2[i][j+1]);
			qr[i][j] = qr[i][j] - f3qt[i][j]*(thp*f2[i][j] + thm*f2[i][j+1]);
		}
	}
		
};

#if defined(ENABLE_ADOLC)
template class Reconstruction<double, adouble>;
template class ReconstructionFirstOrder<double, adouble>;
template class ReconstructionSecondOrder<double, adouble>;
#else
template class Reconstruction<double, double>;
template class ReconstructionFirstOrder<double, double>;
template class ReconstructionSecondOrder<double, double>;
#endif
#endif
