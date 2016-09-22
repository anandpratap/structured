#ifndef _RECONSTRUCTION_H
#define _RECONSTRUCTION_H
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
	virtual void evaluate_chi(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr);

	//! Reconstruct left and right state for all the eta direction faces.
	/*!
	  @param[in] q cell centered variable, includes ghost points. `size: (nic + nhalo) x (njc + nhalo)`
	  @param[out] ql reconstructed left state. `size: nic x nj`
	  @param[out] qr reconstructed right state. `size: nic x nj`
	*/

	virtual void evaluate_eta(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr);
};

/*!
  \brief First order reconstruction

  Simple first order reconstruction where left and right states are taken to be the values of 
  left and right cell.
*/
template<class Tx, class Tad>
class ReconstructionFirstOrder: public Reconstruction<Tx, Tad>{
	size_t ni, nj;
	size_t nic, njc;
public:
	ReconstructionFirstOrder(const size_t val_ni, const size_t val_nj);
	void evaluate_chi(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr);
	void evaluate_eta(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr);
};


/*!
  \brief Second order reconstruction

  Simple second order reconstruction.
*/
template<class Tx, class Tad>
class ReconstructionSecondOrder: public Reconstruction<Tx, Tad>{
	size_t ni, nj;
	size_t nic, njc;
	Tx thm = 2.0/3.0;
	Tx thp = 4.0/3.0;
	Tx eps_chi, eps_eta;
	Array2D<Tad> f2_chi, a1_chi, a2_chi, f3qt_chi;
	Array2D<Tad> f2_eta, a1_eta, a2_eta, f3qt_eta;
public:
	ReconstructionSecondOrder(const size_t val_ni, const size_t val_nj);
	void evaluate_chi(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr);
	void evaluate_eta(const Array2D<const Tad>& q, Array2D<Tad>& ql, Array2D<Tad>& qr);
};
#endif
