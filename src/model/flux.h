#ifndef _FLUX_H
#define _FLUX_H

template <class T>
void roeflux(const T nx, const T ny,
			 const T rlft, const T ulft, const T vlft, const T plft,
			 const T rrht, const T urht, const T vrht, const T prht,
			 T *f){
	constexpr auto gm1 = GAMMA - 1.0;
	constexpr auto ogm1 = 1/gm1;

	const T rlfti = 1/rlft;
	const T rulft = rlft*ulft;
	const T rvlft = rlft*vlft;
	const T uvl = 0.5f*(ulft*ulft + vlft*vlft);
	const T elft = plft*ogm1 + rlft*uvl;
	const T hlft = (elft + plft)*rlfti;


	const T rrhti = 1/rrht;
	const T rurht = rrht*urht;
	const T rvrht = rrht*vrht;
	const T uvr = 0.5f*(urht*urht + vrht*vrht);
	const T erht = prht*ogm1 + rrht*uvr;
	const T hrht = (erht + prht)*rrhti;
  
	const T rat = sqrt(rrht*rlfti);
	const T rati = 1/(rat+1);
	const T rav = rat*rlft;
	const T uav = (rat*urht + ulft)*rati;
	const T vav = (rat*vrht + vlft)*rati;
	const T hav = (rat*hrht + hlft)*rati;
	const T uv = 0.5f*(uav*uav + vav*vav);
	const T cav = sqrt(gm1*(hav - uv));

	T aq1 = rrht - rlft;
	T aq2 = urht - ulft;
	T aq3 = vrht - vlft;
	T aq4 = prht - plft;
  
	const T dr = sqrt(nx*nx + ny*ny);
	const T r1 = nx/dr;
	const T r2 = ny/dr;


	const T uu = r1*uav + r2*vav;
	const T c2 = cav*cav;
	const T c2i = 1/c2;
	const T auu = fabs(uu);
	const T aupc = fabs(uu + cav);
	const T aumc = fabs(uu - cav);

	const T uulft = r1*ulft + r2*vlft;
	const T uurht = r1*urht + r2*vrht;
	const T rcav = rav*cav;
	const T aquu = uurht - uulft;
	const T c2ih = 0.5f*c2i;
	const T ruuav = auu*rav;
  
	const T b1 = auu*(aq1 - c2i*aq4);
	const T b2 = c2ih*aupc*(aq4 + rcav*aquu);
	const T b3 = c2ih*aumc*(aq4 - rcav*aquu);
	const T b4 = b1 + b2 + b3;
	const T b5 = cav*(b2 - b3);
	const T b6 = ruuav*(aq2 - r1*aquu);
	const T b7 = ruuav*(aq3 - r2*aquu);

	aq1 = b4;
	aq2 = uav*b4 + r1*b5 + b6;
	aq3 = vav*b4 + r2*b5 + b7;
	aq4 = hav*b4 + uu*b5 + uav*b6 + vav*b7 - c2*b1*ogm1;

	const T aj = 0.5f*dr;
	const T plar = plft + prht;
	const T eplft = elft + plft;
	const T eprht = erht + prht;
	f[0] = aj*(rlft*uulft + rrht*uurht - aq1);
	f[1] = aj*(rulft*uulft + rurht*uurht + r1*plar - aq2);
	f[2] = aj*(rvlft*uulft + rvrht*uurht + r2*plar - aq3);
	f[3] = aj*(eplft*uulft + eprht*uurht - aq4);
}


#endif
