template <class T>
void roeflux(const T nx, const T ny,
			 const T rlft, const T ulft, const T vlft, const T plft,
			 const T rrht, const T urht, const T vrht, const T prht,
			 T *f){
	T gm1 = GAMMA - 1.0;
	T ogm1 = 1/gm1;

	T rlfti = 1/rlft;
	T rulft = rlft*ulft;
	T rvlft = rlft*vlft;
	T uvl = 0.5f*(ulft*ulft + vlft*vlft);
	T elft = plft*ogm1 + rlft*uvl;
	T hlft = (elft + plft)*rlfti;


	T rrhti = 1/rrht;
	T rurht = rrht*urht;
	T rvrht = rrht*vrht;
	T uvr = 0.5f*(urht*urht + vrht*vrht);
	T erht = prht*ogm1 + rrht*uvr;
	T hrht = (erht + prht)*rrhti;
  
	T rat = sqrt(rrht*rlfti);
	T rati = 1/(rat+1);
	T rav = rat*rlft;
	T uav = (rat*urht + ulft)*rati;
	T vav = (rat*vrht + vlft)*rati;
	T hav = (rat*hrht + hlft)*rati;
	T uv = 0.5f*(uav*uav + vav*vav);
	T cav = sqrt(gm1*(hav - uv));

	T aq1 = rrht - rlft;
	T aq2 = urht - ulft;
	T aq3 = vrht - vlft;
	T aq4 = prht - plft;
  
	T dr = sqrt(nx*nx + ny*ny);
	T r1 = nx/dr;
	T r2 = ny/dr;


	T uu = r1*uav + r2*vav;
	T c2 = cav*cav;
	T c2i = 1/c2;
	T auu = fabs(uu);
	T aupc = fabs(uu + cav);
	T aumc = fabs(uu - cav);

	T uulft = r1*ulft + r2*vlft;
	T uurht = r1*urht + r2*vrht;
	T rcav = rav*cav;
	T aquu = uurht - uulft;
	T c2ih = 0.5f*c2i;
	T ruuav = auu*rav;
  
	T b1, b2, b3, b4, b5, b6, b7;
	b1 = auu*(aq1 - c2i*aq4);
	b2 = c2ih*aupc*(aq4 + rcav*aquu);
	b3 = c2ih*aumc*(aq4 - rcav*aquu);
	b4 = b1 + b2 + b3;
	b5 = cav*(b2 - b3);
	b6 = ruuav*(aq2 - r1*aquu);
	b7 = ruuav*(aq3 - r2*aquu);

	aq1 = b4;
	aq2 = uav*b4 + r1*b5 + b6;
	aq3 = vav*b4 + r2*b5 + b7;
	aq4 = hav*b4 + uu*b5 + uav*b6 + vav*b7 - c2*b1*ogm1;

	T aj = 0.5f*dr;
	T plar = plft + prht;
	T eplft = elft + plft;
	T eprht = erht + prht;
	f[0] = aj*(rlft*uulft + rrht*uurht - aq1);
	f[1] = aj*(rulft*uulft + rurht*uurht + r1*plar - aq2);
	f[2] = aj*(rvlft*uulft + rvrht*uurht + r2*plar - aq3);
	f[3] = aj*(eplft*uulft + eprht*uurht - aq4);
}
