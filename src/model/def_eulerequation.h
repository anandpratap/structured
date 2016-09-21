#ifndef _DEF_EULEREQUATION_H
#define _DEF_EULEREQUATION_H
#include "def_mesh.h"
#include "def_reconstruction.h"
#include "def_flux.h"
#include "def_fluid.h"
#include "def_bc.h"
#include "def_solution.h"
#include "def_io.h"
/*!
  \brief Container for the model EulerEquation
  
  This is a container.
*/

template<class Tx, class Tad>
class EulerEquation{
public:
	size_t ni, nj, nic, njc, nq;
	Array2D<Tad> rho, u, v, p, T, mu, k;
	Array3D<Tad> grad_u, grad_v, grad_T;
	Array2D<Tad> rholft_chi, ulft_chi, vlft_chi, plft_chi;
	Array2D<Tad> rhorht_chi, urht_chi, vrht_chi, prht_chi;
	Array2D<Tad> rholft_eta, ulft_eta, vlft_eta, plft_eta;
	Array2D<Tad> rhorht_eta, urht_eta, vrht_eta, prht_eta;
	Array3D<Tad> flux_chi, flux_eta;

	Array3D<Tad> grad_u_chi, grad_u_eta;
	Array3D<Tad> grad_v_chi, grad_v_eta;
	Array3D<Tad> grad_T_chi, grad_T_eta;
	Array2D<Tad> u_bar_chi, u_bar_eta;
	Array2D<Tad> v_bar_chi, v_bar_eta;
	Array2D<Tad> T_bar_chi, T_bar_eta;
	Array2D<Tad> mu_bar_chi, mu_bar_eta;
	Array2D<Tad> k_bar_chi, k_bar_eta;
	
	std::shared_ptr<Mesh<Tx, Tad>> mesh;
	std::shared_ptr<Config<Tx>> config;
	std::shared_ptr<Reconstruction<Tx, Tad>> reconstruction;

	std::shared_ptr<Reconstruction<Tx, Tad>> reconstruction_rhs;
	std::shared_ptr<Reconstruction<Tx, Tad>> reconstruction_lhs;
	
	std::unique_ptr<ConvectiveFlux<Tx, Tad>> convective_flux;
	std::unique_ptr<DiffusiveFlux<Tx, Tad>> diffusive_flux;
	std::unique_ptr<BoundaryContainer<Tx, Tad>> boundary_container;

	//! calculate residual
	/*!
	  @param a_q[in] flow variable
	  @param a_rhs[out] residual
	*/
	void calc_residual(const Array3D<const Tad>& a_q, Array3D<Tad>& a_rhs, bool lhs=false);
	void calc_dt(const Tx cfl);
	void initialize();
	void calc_convective_residual(Array3D<Tad>& a_rhs);
	void calc_intermediates(const Array3D<const Tad>& a_q);
	void calc_viscous_residual(Array3D<Tad>& a_rhs);
	void calc_source_residual(const Array3D<const Tad>& a_q, Array3D<Tad>& a_rhs);
	void calc_primvars(const Array3D<const Tad>& a_q);
	void calc_boundary();
	EulerEquation(std::shared_ptr<Mesh<Tx, Tad>> val_mesh, std::shared_ptr<Config<Tx>> val_config);
	~EulerEquation(){};
};


#endif
