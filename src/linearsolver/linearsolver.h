#ifndef _LINEARSOLVER_H
#define _LINEARSOLVER_H
#include "common.h"
#include "def_mesh.h"
#include "def_solution.h"

#if defined(ENABLE_ARMA)
#include "def_ls_arma.h"
#endif

#if defined(ENABLE_EIGEN)
#include "def_ls_eigen.h"
#endif

#if defined(ENABLE_PETSC)
#include "petsc.h"
#include "petscksp.h"
#include "def_ls_petsc.h"
#endif

#endif
