#ifndef _DESIGN_H
#define _DESIGN_H
#include "common.h"
template<class Tx, class Tad>
class DesignParameters{
public:
	Array2D<Tx> beta;
	Array2D<Tad> a_beta;
	size_t n;
	DesignParameters(size_t ni, size_t nj);
	DesignParameters(size_t ni);
	~DesignParameters();
};
#endif
