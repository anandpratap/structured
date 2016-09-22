#ifndef _LINEARSOLVER_H
#define _LINEARSOLVER_H
#include "common.h"
#include "mesh.h"
#include "solution.h"

#if defined(ENABLE_ARMA)
#include "ls_arma.h"
#endif

#if defined(ENABLE_EIGEN)
#include "ls_eigen.h"
#endif

#if defined(ENABLE_PETSC)
#include "ls_petsc.h"
#endif

#endif
