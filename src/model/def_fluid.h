#ifndef _DEF_FLUID_H
#define _DEF_FLUID_H
template<class Tx, class Tad>
class FluidModel{
public:
	Tx R, gamma;
	Tx p_ref, rho_ref, T_ref, mu_ref;
	Tx cp, pr;
	
	FluidModel(Tx val_p_ref, Tx val_rho_ref, Tx val_T_ref, Tx val_mu_ref, Tx val_pr = 0.7);
	~FluidModel(){};

	template<class Tq>
	inline Tq get_T_prho(const Tq p, const Tq rho);

	template<class Tq>
	inline Tq get_rho_pT(const Tq p, const Tq T);

	template<class Tq>
	inline Tq get_p_rhoT(const Tq rho, const Tq T);

	template<class Tq>
	inline Tq get_laminar_viscosity(const Tq T);

	template<class Tq>
	inline Tq get_thermal_conductivity(const Tq T);

	template<class Tq>
	void primvars(const Array3D<const Tq>& Q, Array2D<Tq>& rho, Array2D<Tq>& u, Array2D<Tq>& v, Array2D<Tq>& p, Array2D<Tq>& T, const size_t shifti = 0, const size_t shiftj = 0);
};
#endif
