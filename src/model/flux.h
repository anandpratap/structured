#ifndef _FLUX_H
#define _FLUX_H
template <class T, class Tad>
class ConvectiveFlux{
public:
	virtual void evaluate(Array3D<T>& normal,
				  Array2D<Tad>& rlft_a, Array2D<Tad>& ulft_a, Array2D<Tad>& vlft_a, Array2D<Tad>& plft_a,
				  Array2D<Tad>& rrht_a, Array2D<Tad>& urht_a, Array2D<Tad>& vrht_a, Array2D<Tad>& prht_a,
				  Array3D<Tad>& f_a){};
};

template<class T, class Tad>
class RoeFlux: public ConvectiveFlux<T, Tad>{
 public:
	void evaluate(Array3D<T>& normal,
				  Array2D<Tad>& rlft_a, Array2D<Tad>& ulft_a, Array2D<Tad>& vlft_a, Array2D<Tad>& plft_a,
				  Array2D<Tad>& rrht_a, Array2D<Tad>& urht_a, Array2D<Tad>& vrht_a, Array2D<Tad>& prht_a,
				  Array3D<Tad>& f_a){
		constexpr auto gm1 = GAMMA - 1.0;
		constexpr auto ogm1 = 1.0/gm1;
		auto ni = normal.extent(0);
		auto nj = normal.extent(1);

#pragma omp parallel for		
		for(int i=0; i<ni; i++){
			for(int j=0; j<nj; j++){
				const T nx = normal[i][j][0];
				const T ny = normal[i][j][1];
				const Tad rlft = rlft_a[i][j];
				const Tad ulft = ulft_a[i][j];
				const Tad vlft = vlft_a[i][j];
				const Tad plft = plft_a[i][j];

				const Tad rrht = rrht_a[i][j];
				const Tad urht = urht_a[i][j];
				const Tad vrht = vrht_a[i][j];
				const Tad prht = prht_a[i][j];

				const Tad rlfti = 1.0/rlft;
				const Tad rulft = rlft*ulft;
				const Tad rvlft = rlft*vlft;
				const Tad uvl = 0.5*(ulft*ulft + vlft*vlft);
				const Tad elft = plft*ogm1 + rlft*uvl;
				const Tad hlft = (elft + plft)*rlfti;
		
				
				const Tad rrhti = 1.0/rrht;
				const Tad rurht = rrht*urht;
				const Tad rvrht = rrht*vrht;
				const Tad uvr = 0.5*(urht*urht + vrht*vrht);
				const Tad erht = prht*ogm1 + rrht*uvr;
				const Tad hrht = (erht + prht)*rrhti;
				
				const Tad rat = sqrt(rrht*rlfti);
				const Tad rati = 1.0/(rat+1.0);
				const Tad rav = rat*rlft;
				const Tad uav = (rat*urht + ulft)*rati;
				const Tad vav = (rat*vrht + vlft)*rati;
				const Tad hav = (rat*hrht + hlft)*rati;
				const Tad uv = 0.5*(uav*uav + vav*vav);
				const Tad cav = sqrt(gm1*(hav - uv));
				
				Tad aq1 = rrht - rlft;
				Tad aq2 = urht - ulft;
				Tad aq3 = vrht - vlft;
				Tad aq4 = prht - plft;
				
				const Tad dr = sqrt(nx*nx + ny*ny);
				const Tad r1 = nx/dr;
				const Tad r2 = ny/dr;
				
				
				const Tad uu = r1*uav + r2*vav;
				const Tad c2 = cav*cav;
				const Tad c2i = 1.0/c2;
				const Tad auu = fabs(uu);
				const Tad aupc = fabs(uu + cav);
				const Tad aumc = fabs(uu - cav);
				
				const Tad uulft = r1*ulft + r2*vlft;
				const Tad uurht = r1*urht + r2*vrht;
				const Tad rcav = rav*cav;
				const Tad aquu = uurht - uulft;
				const Tad c2ih = 0.5*c2i;
				const Tad ruuav = auu*rav;
				
				const Tad b1 = auu*(aq1 - c2i*aq4);
				const Tad b2 = c2ih*aupc*(aq4 + rcav*aquu);
				const Tad b3 = c2ih*aumc*(aq4 - rcav*aquu);
				const Tad b4 = b1 + b2 + b3;
				const Tad b5 = cav*(b2 - b3);
				const Tad b6 = ruuav*(aq2 - r1*aquu);
				const Tad b7 = ruuav*(aq3 - r2*aquu);
				
				aq1 = b4;
				aq2 = uav*b4 + r1*b5 + b6;
				aq3 = vav*b4 + r2*b5 + b7;
				aq4 = hav*b4 + uu*b5 + uav*b6 + vav*b7 - c2*b1*ogm1;
				
				const Tad aj = 0.5*dr;
				const Tad plar = plft + prht;
				const Tad eplft = elft + plft;
				const Tad eprht = erht + prht;
				f_a[i][j][0] = aj*(rlft*uulft + rrht*uurht - aq1);
				f_a[i][j][1] = aj*(rulft*uulft + rurht*uurht + r1*plar - aq2);
				f_a[i][j][2] = aj*(rvlft*uulft + rvrht*uurht + r2*plar - aq3);
				f_a[i][j][3] = aj*(eplft*uulft + eprht*uurht - aq4);
			}
		}
	};
};


template<class T, class Tad>
class AUSMFlux: public ConvectiveFlux<T, Tad>{
 public:
	void evaluate(Array3D<T>& normal,
				  Array2D<Tad>& rlft_a, Array2D<Tad>& ulft_a, Array2D<Tad>& vlft_a, Array2D<Tad>& plft_a,
				  Array2D<Tad>& rrht_a, Array2D<Tad>& urht_a, Array2D<Tad>& vrht_a, Array2D<Tad>& prht_a,
				  Array3D<Tad>& f_a){
		constexpr auto gm1 = GAMMA - 1.0;
		constexpr auto ogm1 = 1.0/gm1;
		auto ni = normal.extent(0);
		auto nj = normal.extent(1);

#pragma omp parallel for		
		for(int i=0; i<ni; i++){
			for(int j=0; j<nj; j++){
				const T nx = normal[i][j][0];
				const T ny = normal[i][j][1];
				const Tad rlft = rlft_a[i][j];
				const Tad ulft = ulft_a[i][j];
				const Tad vlft = vlft_a[i][j];
				const Tad plft = plft_a[i][j];

				const Tad rrht = rrht_a[i][j];
				const Tad urht = urht_a[i][j];
				const Tad vrht = vrht_a[i][j];
				const Tad prht = prht_a[i][j];
				
				const T ds = std::sqrt(nx*nx + ny*ny);
				const Tad ulft_normal = (ulft*nx + vlft*ny)/ds;
				const Tad urht_normal = (urht*nx + vrht*ny)/ds;

				const Tad alft = std::sqrt(GAMMA*plft/rlft);
				const Tad arht = std::sqrt(GAMMA*prht/rrht);

				const Tad machlft = ulft_normal/alft;
				const Tad machrht = urht_normal/arht;
				

				const Tad rlfti = 1.0/rlft;
				const Tad uvl = 0.5*(ulft*ulft + vlft*vlft);
				const Tad elft = plft*ogm1 + rlft*uvl;
				const Tad hlft = (elft + plft)*rlfti;
		
				
				const Tad rrhti = 1.0/rrht;
				const Tad uvr = 0.5*(urht*urht + vrht*vrht);
				const Tad erht = prht*ogm1 + rrht*uvr;
				const Tad hrht = (erht + prht)*rrhti;
				
				auto mach_p = [](Tad M){return std::abs(M) <= 1.0 ? 0.25*(M+1.0)*(M+1.0): 0.5*(M + fabs(M));};
				auto mach_m = [](Tad M){return std::abs(M) <= 1.0 ? -0.25*(M-1.0)*(M-1.0): 0.5*(M - fabs(M));};

				const Tad mach_half = mach_p(machlft) + mach_m(machrht);

				auto pres_p = [](Tad M, Tad p){return std::abs(M) <= 1.0 ? 0.25*p*(M+1.0)*(M+1.0)*(2.0-M): 0.5*p*(M+fabs(M))/M;};
				auto pres_m = [](Tad M, Tad p){return std::abs(M) <= 1.0 ? 0.25*p*(M-1.0)*(M-1.0)*(2.0+M): 0.5*p*(M-fabs(M))/M;};
				const Tad p_half = pres_p(machlft, plft) + pres_m(machrht, prht);

				if(mach_half >= 0.0){
					f_a[i][j][0] = rlft*alft*mach_half*ds;
					f_a[i][j][1] = rlft*alft*ulft*mach_half*ds + p_half*nx;
					f_a[i][j][2] = rlft*alft*vlft*mach_half*ds + p_half*ny;
					f_a[i][j][3] = rlft*alft*hlft*mach_half*ds;
				}
				else{
					f_a[i][j][0] = rrht*arht*mach_half*ds;
					f_a[i][j][1] = rrht*arht*urht*mach_half*ds + p_half*nx;
					f_a[i][j][2] = rrht*arht*vrht*mach_half*ds + p_half*ny;
					f_a[i][j][3] = rrht*arht*hrht*mach_half*ds;
				}


				
			}
		}
	};
};

#endif
