#ifndef _MESH_H
#define _MESH_H
#include "common.h"
#include "config.h"
#include <memory>
#include "fluid.h"
template<class Tx, class Tad>
class EulerEquation;
template<class Tx, class Tad>
class Solution;
template<class Tx, class Tad>
class IOManager;

template<class Tx, class Tad>
class LinearSolverArma;

template<class Tx, class Tad>
class LinearSolverEigen;

template<class Tx, class Tad>
class LinearSolverPetsc;

template <class Tx, class Tad>
class Mesh: public std::enable_shared_from_this<Mesh<Tx, Tad>>{
 public:
	std::string label;
	std::shared_ptr<Config<Tx>> config;
	size_t ni, nj;
	size_t nic, njc;
	size_t j1, nb;
	Array2D<Tx> xv, yv;
	Array2D<Tx> xc, yc;

	Array3D<Tx> normal_eta;
	Array3D<Tx> normal_chi;

	Array2D<Tx> ds_eta;
	Array2D<Tx> ds_chi;

	Array2D<Tx> volume;
	Array2D<Tx> chi_x, chi_y, eta_x, eta_y;
	Array2D<Tx> x_chi, y_chi, x_eta, y_eta;
	
	std::shared_ptr<Solution<Tx,Tad>> solution;
	std::shared_ptr<IOManager<Tx, Tad>> iomanager;
	std::shared_ptr<EulerEquation<Tx, Tad>> equation;
	std::shared_ptr<FluidModel<Tx, Tad>> fluid_model;

#if defined(ENABLE_ARMA)
	std::shared_ptr<LinearSolverArma<Tx, Tad>> linearsolver;
#endif
#if defined(ENABLE_EIGEN)
	std::shared_ptr<LinearSolverEigen<Tx, Tad>> linearsolver;
#endif
#if defined(ENABLE_PETSC)
	std::shared_ptr<LinearSolverPetsc<Tx, Tad>> linearsolver;
#endif

 public:
	Mesh(std::shared_ptr<Config<Tx>> config);

	Mesh(std::shared_ptr<Mesh<Tx, Tad>> mesh, const size_t nskipi=0, const size_t nskipj=0, const size_t refine=0);
	~Mesh();

	void calc_metrics();

	template<class Tq>
		void calc_gradient(const Array2D<const Tq>& q, Array3D<Tq> &grad_q, size_t skipi=0, size_t skipj=0);

	template<class Tq>
		void calc_gradient(const Array2D<const Tq>& q, Array3D<Tq>& grad_chi, Array3D<Tq>& grad_eta);

	template<class Tq>
		void calc_face(const Array2D<const Tq>& q, Array2D<Tq>& q_chi, Array2D<Tq>& q_eta);	
	void setup();

	void simple_loader(std::string filename);
	void plot3d_loader(std::string filename);
};

#endif
