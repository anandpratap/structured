#ifndef _DESIGN_CPP
#define _DESIGN_CPP
#include "design.h"
template<class Tx, class Tad>
DesignParameters<Tx, Tad>::DesignParameters(size_t ni, size_t nj){
	n = ni*nj;
	beta = Array2D<Tx>(ni, nj);
	a_beta = Array2D<Tad>(ni, nj);
	beta.fill(1.0+1e-2);
	a_beta.fill(1.0+1e-2);
}

template<class Tx, class Tad>
DesignParameters<Tx, Tad>::DesignParameters(size_t ni){
	n = ni;
	beta = Array2D<Tx>(ni, 1);
	a_beta = Array2D<Tad>(ni, 1);
	beta.fill(1.0);
	a_beta.fill(1.0);
}

template<class Tx, class Tad>
DesignParameters<Tx, Tad>::~DesignParameters(){};

#if defined(ENABLE_ADOLC)
template class DesignParameters<double, adouble>;
#else
#endif

#endif
