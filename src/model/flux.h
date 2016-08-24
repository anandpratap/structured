#ifndef _FLUX_H
#define _FLUX_H
template <class T, class Tad>
	void roeflux(const T nx, const T ny,
				 const Tad rlft, const Tad ulft, const Tad vlft, const Tad plft,
				 const Tad rrht, const Tad urht, const Tad vrht, const Tad prht,
				 Tad *f){
	constexpr auto gm1 = GAMMA - 1.0;
	constexpr auto ogm1 = 1/gm1;

	const Tad rlfti = 1/rlft;
	const Tad rulft = rlft*ulft;
	const Tad rvlft = rlft*vlft;
	const Tad uvl = 0.5f*(ulft*ulft + vlft*vlft);
	const Tad elft = plft*ogm1 + rlft*uvl;
	const Tad hlft = (elft + plft)*rlfti;


	const Tad rrhti = 1/rrht;
	const Tad rurht = rrht*urht;
	const Tad rvrht = rrht*vrht;
	const Tad uvr = 0.5f*(urht*urht + vrht*vrht);
	const Tad erht = prht*ogm1 + rrht*uvr;
	const Tad hrht = (erht + prht)*rrhti;
  
	const Tad rat = sqrt(rrht*rlfti);
	const Tad rati = 1/(rat+1);
	const Tad rav = rat*rlft;
	const Tad uav = (rat*urht + ulft)*rati;
	const Tad vav = (rat*vrht + vlft)*rati;
	const Tad hav = (rat*hrht + hlft)*rati;
	const Tad uv = 0.5f*(uav*uav + vav*vav);
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
	const Tad c2i = 1/c2;
	const Tad auu = fabs(uu);
	const Tad aupc = fabs(uu + cav);
	const Tad aumc = fabs(uu - cav);

	const Tad uulft = r1*ulft + r2*vlft;
	const Tad uurht = r1*urht + r2*vrht;
	const Tad rcav = rav*cav;
	const Tad aquu = uurht - uulft;
	const Tad c2ih = 0.5f*c2i;
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

	const Tad aj = 0.5f*dr;
	const Tad plar = plft + prht;
	const Tad eplft = elft + plft;
	const Tad eprht = erht + prht;
	f[0] = aj*(rlft*uulft + rrht*uurht - aq1);
	f[1] = aj*(rulft*uulft + rurht*uurht + r1*plar - aq2);
	f[2] = aj*(rvlft*uulft + rvrht*uurht + r2*plar - aq3);
	f[3] = aj*(eplft*uulft + eprht*uurht - aq4);
}


#endif
