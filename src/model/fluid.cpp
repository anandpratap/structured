#ifndef _FLUID_H
#define _FLUID_H
#include "def_fluid.h"
template<class Tx, class Tad>
FluidModel<Tx, Tad>::FluidModel(Tx val_p_ref, Tx val_rho_ref, Tx val_T_ref, Tx val_mu_ref, Tx val_pr){
	p_ref = val_p_ref;
	rho_ref = val_rho_ref;
	T_ref = val_T_ref;
	mu_ref = val_mu_ref;
	pr = val_pr;
		
	R = p_ref/rho_ref/T_ref;
	gamma = 1.4;
	cp = gamma*R/(gamma-1.0);
};

template<class Tx, class Tad>
template<class Tq>
Tq FluidModel<Tx, Tad>::get_T_prho(const Tq p, const Tq rho){
	return p/rho/R;
};

template<class Tx, class Tad>
template<class Tq>
Tq FluidModel<Tx, Tad>::get_rho_pT(const Tq p, const Tq T){
	return p/T/R;
};


template<class Tx, class Tad>
template<class Tq>
Tq FluidModel<Tx, Tad>::get_p_rhoT(const Tq rho, const Tq T){
	return rho*R*T;
};

template<class Tx, class Tad>
template<class Tq>
Tq FluidModel<Tx, Tad>::get_laminar_viscosity(const Tq T){
	return mu_ref*pow(T/T_ref, 2.0/3.0);
};

template<class Tx, class Tad>
template<class Tq>
Tq FluidModel<Tx, Tad>::get_thermal_conductivity(const Tq T){
	return get_laminar_viscosity(T)*cp/pr;
};

template<class Tx, class Tad>
template<class Tq>
void FluidModel<Tx, Tad>::primvars(const Array3D<const Tq>& Q, Array2D<Tq>& rho, Array2D<Tq>& u, Array2D<Tq>& v, Array2D<Tq>& p, Array2D<Tq>& T, const size_t shifti, const size_t shiftj){
	auto nic = Q.extent(0);
	auto njc = Q.extent(1);
#pragma omp parallel for
	for(size_t i=0; i<nic; i++){
		for(size_t j=0; j<njc; j++){
			Tq tmp_rho, tmp_u, tmp_v;
			tmp_rho = Q[i][j][0];
			tmp_u = Q[i][j][1]/tmp_rho;
			tmp_v = Q[i][j][2]/tmp_rho;
			rho[i+shifti][j+shiftj] = tmp_rho;
			u[i+shifti][j+shiftj] = tmp_u;
			v[i+shifti][j+shiftj] = tmp_v;
			p[i+shifti][j+shiftj] = (Q[i][j][3] - 0.5*tmp_rho*(tmp_u*tmp_u + tmp_v*tmp_v))*(gamma-1.0);
			T[i+shifti][j+shiftj] = get_T_prho(p[i+shifti][j+shiftj], rho[i+shifti][j+shiftj]);
		}
	}
}
#if defined(ENABLE_ADOLC)
template class FluidModel<double,adouble>;
template adouble FluidModel<double,adouble>::get_p_rhoT<adouble>(const adouble, const adouble);
template adouble FluidModel<double,adouble>::get_T_prho<adouble>(const adouble, const adouble);
template adouble FluidModel<double,adouble>::get_rho_pT<adouble>(const adouble, const adouble);
template adouble FluidModel<double,adouble>::get_laminar_viscosity<adouble>(const adouble);
template adouble FluidModel<double,adouble>::get_thermal_conductivity<adouble>(const adouble);
template void FluidModel<double,adouble>::primvars<adouble>(const Array3D<const adouble>& Q, Array2D<adouble>& rho, Array2D<adouble>& u, Array2D<adouble>& v, Array2D<adouble>& p, Array2D<adouble>& T, const size_t shifti, const size_t shiftj);

template double FluidModel<double,adouble>::get_p_rhoT<double>(const double, const double);
template double FluidModel<double,adouble>::get_T_prho<double>(const double, const double);
template double FluidModel<double,adouble>::get_rho_pT<double>(const double, const double);
template double FluidModel<double,adouble>::get_laminar_viscosity<double>(const double);
template double FluidModel<double,adouble>::get_thermal_conductivity<double>(const double);
template void FluidModel<double,adouble>::primvars<double>(const Array3D<const double>& Q, Array2D<double>& rho, Array2D<double>& u, Array2D<double>& v, Array2D<double>& p, Array2D<double>& T, const size_t shifti, const size_t shiftj);
#else
template class FluidModel<double,double>;
template class FluidModel<float,float>;
#endif

#endif
