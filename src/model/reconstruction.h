#ifndef __RECONSTRUCTION__H
#define __RECONSTRUCTION__H
#include "common.h"

/*!
  \brief Base class for classes used to reconstruct/interpolate left and right state for a riemann solver.
*/
template<class Tx, class Tad>
class Reconstruction{
public:
	//! Reconstruct left and right state for all the chi direction faces.
	/*!
	  @param[in] q cell centered variable, includes ghost points. `size: (nic + nhalo) x (njc + nhalo)`
	  @param[out] ql reconstructed left state. `size: ni x njc`
	  @param[out] qr reconstructed right state. `size: ni x njc`
	*/
	virtual void evaluate_chi(Array2D<Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
	};

	//! Reconstruct left and right state for all the eta direction faces.
	/*!
	  @param[in] q cell centered variable, includes ghost points. `size: (nic + nhalo) x (njc + nhalo)`
	  @param[out] ql reconstructed left state. `size: nic x nj`
	  @param[out] qr reconstructed right state. `size: nic x nj`
	*/

	virtual void evaluate_eta(Array2D<Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
	};
};

/*!
  \brief First order reconstruction

  Simple first order reconstruction where left and right states are taken to be the values of 
  left and right cell.
 */
template<class Tx, class Tad>
class ReconstructionFirstOrder: public Reconstruction<Tx, Tad>{
	uint ni, nj;
	uint nic, njc;
public:
 ReconstructionFirstOrder(const uint val_ni, const uint val_nj): Reconstruction<Tx, Tad>(){
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
	void evaluate_chi(Array2D<Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
#pragma omp parallel for
		for(uint i=0; i<ni; i++){
			for(uint j=0; j<njc; j++){
				ql[i][j] = q[i][j+1];
				qr[i][j] = q[i+1][j+1];
			}
		}
	};
	
	void evaluate_eta(Array2D<Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
#pragma omp parallel for
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<nj; j++){
				ql[i][j] = q[i+1][j];
				qr[i][j] = q[i+1][j+1];
			}
		}
	};
};


/*!
  \brief Second order reconstruction

  Simple second order reconstruction.
 */
template<class Tx, class Tad>
class ReconstructionSecondOrder: public Reconstruction<Tx, Tad>{
	uint ni, nj;
	uint nic, njc;
	Tx thm = 2.0/3.0;
	Tx thp = 4.0/3.0;
	Tx eps_chi, eps_eta;
	Array2D<Tad> f2_chi, a1_chi, a2_chi, f3qt_chi;
	Array2D<Tad> f2_eta, a1_eta, a2_eta, f3qt_eta;
public:
 ReconstructionSecondOrder(const uint val_ni, const uint val_nj): Reconstruction<Tx, Tad>(){
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
	void evaluate_chi(Array2D<Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
		auto f2 = f2_chi;
		auto a1 = a1_chi;
		auto a2 = a2_chi;
		auto f3qt = f3qt_chi;
		auto eps = eps_chi;
#pragma omp parallel  for
		for(uint i=0; i<ni; i++){
			for(uint j=0; j<njc; j++){
				ql[i][j] = q[i][j+1];
				qr[i][j] = q[i+1][j+1];
			}
		}
#pragma omp parallel  for
		for(uint i=0; i<ni; i++){
			for(uint j=0; j<njc; j++){
				f2[i][j] = q[i+1][j+1] - q[i][j+1];
			}
		}
		
#pragma omp parallel  for
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				a1[i][j] = 3.0*f2[i+1][j]*f2[i][j];
				a2[i][j] = 2.0*(f2[i+1][j] - f2[i][j])*(f2[i+1][j] - f2[i][j]) + a1[i][j];
				f3qt[i][j] = 0.25*(a1[i][j] + eps)/(a2[i][j] + eps);
			}
		}
		
#pragma omp parallel  for
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				ql[i+1][j] = ql[i+1][j] + f3qt[i][j]*(thm*f2[i][j] + thp*f2[i+1][j]);
				qr[i][j] = qr[i][j] - f3qt[i][j]*(thp*f2[i][j] + thm*f2[i+1][j]);
			}
		}
	};
	
	void evaluate_eta(Array2D<Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr){
		auto f2 = f2_eta;
		auto a1 = a1_eta;
		auto a2 = a2_eta;
		auto f3qt = f3qt_eta;
		auto eps = eps_eta;
#pragma omp parallel  for
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<nj; j++){
				ql[i][j] = q[i+1][j];
				qr[i][j] = q[i+1][j+1];
			}
		}
#pragma omp parallel  for
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<nj; j++){
				f2[i][j] = q[i+1][j+1] - q[i+1][j];
			}
		}
		
#pragma omp parallel  for
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				a1[i][j] = 3.0*f2[i][j+1]*f2[i][j];
				a2[i][j] = 2.0*(f2[i][j+1] - f2[i][j])*(f2[i][j+1] - f2[i][j]) + a1[i][j];
				f3qt[i][j] = 0.25*(a1[i][j] + eps)/(a2[i][j] + eps);
			}
		}
		
#pragma omp parallel  for
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
				ql[i][j+1] = ql[i][j+1] + f3qt[i][j]*(thm*f2[i][j] + thp*f2[i][j+1]);
				qr[i][j] = qr[i][j] - f3qt[i][j]*(thp*f2[i][j] + thm*f2[i][j+1]);
			}
		}
		
	};
};
#endif
