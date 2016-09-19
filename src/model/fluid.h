#ifndef _GAS_H
#define _GAS_H
template<class Tx, class Tad>
class FluidModel{
public:
	Tx R, gamma;
	Tx p_ref, rho_ref, T_ref, mu_ref;
	Tx cp, pr;
	
	FluidModel(Tx val_p_ref, Tx val_rho_ref, Tx val_T_ref, Tx val_mu_ref, Tx val_pr = 0.7){
		p_ref = val_p_ref;
		rho_ref = val_rho_ref;
		T_ref = val_T_ref;
		mu_ref = val_mu_ref;
		pr = val_pr;
		
		R = p_ref/rho_ref/T_ref;
		gamma = 1.4;
		cp = gamma*R/(gamma-1.0);
	};
	~FluidModel(){};

	template<class Tq>
	inline Tq get_T_prho(Tq p, Tq rho){
		return p/rho/R;
	};

	template<class Tq>
		inline Tq get_rho_pT(Tq p, Tq T){
		return p/T/R;
	};

	template<class Tq>
	inline Tq get_p_rhoT(Tq rho, Tq T){
		return rho*R*T;
	};

	template<class Tq>
	inline Tq get_laminar_viscosity(Tq T){
		return mu_ref*pow(T/T_ref, 2.0/3.0);
	};

	template<class Tq>
	inline Tq get_thermal_conductivity(Tq T){
		return get_laminar_viscosity(T)*cp/pr;
	};

	template<class Tq>
	void primvars(Array3D<Tq>& Q, Array2D<Tq>& rho, Array2D<Tq>& u, Array2D<Tq>& v, Array2D<Tq>& p, Array2D<Tq>& T, const uint shifti = 0, const uint shiftj = 0){
		uint nic = Q.extent(0);
		uint njc = Q.extent(1);
#pragma omp parallel for
		for(uint i=0; i<nic; i++){
			for(uint j=0; j<njc; j++){
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
};
#endif
