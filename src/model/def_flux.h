#ifndef _DEF_FLUX_H
#define _DEF_FLUX_H
#include "common.h"

template <class Tx, class Tad>
class DiffusiveFlux{
public:
	virtual void evaluate(const Array3D<const Tx>& val_normal,
						  const Array3D<const Tad>& val_grad_u, const Array3D<const Tad>& val_grad_v, const Array3D<const Tad>& val_grad_T,
						  const Array2D<const Tad>& val_ubar, const Array2D<const Tad>& val_vbar,
						  const Array2D<const Tad>& val_mubar, const Array2D<const Tad>& val_kbar,
						  Array3D<Tad>& val_flux);
};

template<class Tx, class Tad>
class DiffusiveFluxGreenGauss: public DiffusiveFlux<Tx, Tad>{
public:
	virtual void evaluate(const Array3D<const Tx>& val_normal,
						  const Array3D<const Tad>& val_grad_u, const Array3D<const Tad>& val_grad_v, const Array3D<const Tad>& val_grad_T,
						  const Array2D<const Tad>& val_ubar, const Array2D<const Tad>& val_vbar,
						  const Array2D<const Tad>& val_mubar, const Array2D<const Tad>& val_kbar,
						  Array3D<Tad>& val_flux);
};

template <class Tx, class Tad>
class ConvectiveFlux{
public:
	virtual void evaluate(const Array3D<const Tx>& normal,
						  const Array2D<const Tad>& rlft_a, const Array2D<const Tad>& ulft_a, const Array2D<const Tad>& vlft_a, const Array2D<const Tad>& plft_a,
						  const Array2D<const Tad>& rrht_a, const Array2D<const Tad>& urht_a, const Array2D<const Tad>& vrht_a, const Array2D<const Tad>& prht_a,
						  Array3D<Tad>& f_a){};
};

template<class Tx, class Tad>
class ConvectiveFluxRoe: public ConvectiveFlux<Tx, Tad>{
public:
	void evaluate(const Array3D<const Tx>& normal,
				  const Array2D<const Tad>& rlft_a, const Array2D<const Tad>& ulft_a, const Array2D<const Tad>& vlft_a, const Array2D<const Tad>& plft_a,
				  const Array2D<const Tad>& rrht_a, const Array2D<const Tad>& urht_a, const Array2D<const Tad>& vrht_a, const Array2D<const Tad>& prht_a,
				  Array3D<Tad>& f_a);
};


template<class Tx, class Tad>
class ConvectiveFluxAUSM: public ConvectiveFlux<Tx, Tad>{
public:
	Tad mach_p(const Tad M);
	Tad mach_m(const Tad M);
	Tad pres_p(const Tad M, const Tad p);
	Tad pres_m(const Tad M, const Tad p);
			
	void evaluate(const Array3D<const Tx>& normal,
				  const Array2D<const Tad>& rlft_a, const Array2D<const Tad>& ulft_a, const Array2D<const Tad>& vlft_a, const Array2D<const Tad>& plft_a,
				  const Array2D<const Tad>& rrht_a, const Array2D<const Tad>& urht_a, const Array2D<const Tad>& vrht_a, const Array2D<const Tad>& prht_a,
				  Array3D<Tad>& f_a);
};

#endif
