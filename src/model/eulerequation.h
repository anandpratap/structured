#ifndef _EULEREQUATION_H
#define _EULEREQUATION_H
#include "flux.h"
#include "reconstruction.h"
#include "bc.h"

template<class Tx, class Tad>
class EulerEquation{
public:
	uint ni, nj, nic, njc, nq;
	Array2D<Tad> rho, u, v, p, T;
	Array3D<Tad> grad_u, grad_v, grad_T;
	Array2D<Tad> rholft_xi, ulft_xi, vlft_xi, plft_xi;
	Array2D<Tad> rhorht_xi, urht_xi, vrht_xi, prht_xi;
	Array2D<Tad> rholft_eta, ulft_eta, vlft_eta, plft_eta;
	Array2D<Tad> rhorht_eta, urht_eta, vrht_eta, prht_eta;
	Array3D<Tad> flux_xi, flux_eta;

	std::shared_ptr<Mesh<Tx>> mesh;
	std::shared_ptr<Config<Tx>> config;
	std::unique_ptr<Reconstruction<Tx, Tad>> reconstruction;
	std::unique_ptr<ConvectiveFlux<Tx, Tad>> convective_flux;
	

	std::vector<BoundaryCondition<Tx, Tad>*> boundaryconditions;
	//std::unique_ptr<TransportEquation<Tx, Tad>> transport;	
	void calc_residual(Array3D<Tad>& a_q, Array3D<Tad>& a_rhs);
	void calc_convective_residual(Array3D<Tad>& a_rhs);
	void calc_viscous_residual(Array3D<Tad>& a_rhs){
		auto mu = config->freestream->mu_inf;
		auto p_inf = config->freestream->p_inf;
		auto pr_inf = config->freestream->pr_inf;
		auto rho_inf = config->freestream->rho_inf;
		auto T_inf = config->freestream->T_inf;
		auto Rc = p_inf/rho_inf/T_inf;
		auto k = mu*GAMMA*Rc/(GAMMA-1.0)/pr_inf;
		mesh->calc_gradient(u, grad_u, 1, 1);
		mesh->calc_gradient(v, grad_v, 1, 1);
		mesh->calc_gradient(T, grad_T, 1, 1);
		for(int i=0; i<nic; i++){
			for(int j=0; j<njc; j++){
				mesh->solution->q_aux[i][j][0] = grad_u[i][j][0].value();
				mesh->solution->q_aux[i][j][1] = grad_u[i][j][1].value();
				mesh->solution->q_aux[i][j][2] = grad_v[i][j][0].value();
				mesh->solution->q_aux[i][j][3] = grad_v[i][j][1].value();
			}
		}
		
		// xi
		
		for(int i=0; i<ni; i++){
			for(int j=0; j<njc; j++){
				Tad dudx, dudy, dvdx, dvdy, ubar, vbar, dTdx, dTdy;
				auto nx = mesh->normal_chi[i][j][0];
					auto ny = mesh->normal_chi[i][j][1];
					if(i == 0){
											/* dvdy = grad_v[i][j][1]; */
					const uint b = 1;

					const Tad uleft = u[i+b-1][j+b];
					const Tad uright = u[i+b][j+b];
					const Tad utop = 0.25*(uleft+uright+u[i+b-1][j+b+1]+u[i+b][j+b+1]);
					const Tad ubottom = 0.25*(uleft+uright+u[i+b-1][j+b-1]+u[i+b][j+b-1]);

					const Tad vleft = v[i+b-1][j+b];
					const Tad vright = v[i+b][j+b];
					const Tad vtop = 0.25*(vleft+vright+v[i+b-1][j+b+1]+v[i+b][j+b+1]);
					const Tad vbottom = 0.25*(vleft+vright+v[i+b-1][j+b-1]+v[i+b][j+b-1]);


					
					const Tad Tleft = T[i+b-1][j+b];
					const Tad Tright = T[i+b][j+b];
					const Tad Ttop = 0.25*(Tleft+Tright+T[i+b-1][j+b+1]+T[i+b][j+b+1]);
					const Tad Tbottom = 0.25*(Tleft+Tright+T[i+b-1][j+b-1]+T[i+b][j+b-1]);
					ubar = 0.25*(uleft+uright+utop+ubottom);
					vbar = 0.25*(vleft+vright+vtop+vbottom);
					

						
						double normal_top[2];
						double normal_bottom[2];
						double normal_left[2];
						double normal_right[2];
						// for(int k=0; k<2; k++){
						// 	normal_top[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j+1][k]);
						// 	normal_bottom[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j-1][k]);
						// 	normal_left[k] = 0.5*(mesh->normal_chi[i][j] + mesh->normal_chi[i][j-1]);
						// 	normal_right[k] = 0.5*(mesh->normal_chi[i+1][j] + mesh->normal_chi[i+1][j-1]);
						// }
						double volume = mesh->volume[i][j];
						for(int k=0; k<2; k++){
							normal_top[k] = mesh->normal_eta[i][j+1][k];
							normal_bottom[k] = mesh->normal_eta[i][j][k];
							normal_right[k] = 0.5*(mesh->normal_chi[i][j][k] + mesh->normal_chi[i+1][j][k]);
							normal_left[k] = mesh->normal_chi[i][j][k];
						}

						
						dudx = (normal_top[0]*utop - normal_bottom[0]*ubottom + normal_right[0]*uright - normal_left[0]*uleft)/volume;
						dudy = (normal_top[1]*utop - normal_bottom[1]*ubottom + normal_right[1]*uright - normal_left[1]*uleft)/volume;
						
						dvdx = (normal_top[0]*vtop - normal_bottom[0]*vbottom + normal_right[0]*vright - normal_left[0]*vleft)/volume;
						dvdy = (normal_top[1]*vtop - normal_bottom[1]*vbottom + normal_right[1]*vright - normal_left[1]*vleft)/volume;
						
						dTdx = (normal_top[0]*Ttop - normal_bottom[0]*Tbottom + normal_right[0]*Tright - normal_left[0]*Tleft)/volume;
						dTdy = (normal_top[1]*Ttop - normal_bottom[1]*Tbottom + normal_right[1]*Tright - normal_left[1]*Tleft)/volume;

				}

				else if(i == ni-1){
					/* dvdy = grad_v[i][j][1]; */
					const uint b = 1;

					const Tad uleft = u[i+b-1][j+b];
					const Tad uright = u[i+b][j+b];
					const Tad utop = 0.25*(uleft+uright+u[i+b-1][j+b+1]+u[i+b][j+b+1]);
					const Tad ubottom = 0.25*(uleft+uright+u[i+b-1][j+b-1]+u[i+b][j+b-1]);

					const Tad vleft = v[i+b-1][j+b];
					const Tad vright = v[i+b][j+b];
					const Tad vtop = 0.25*(vleft+vright+v[i+b-1][j+b+1]+v[i+b][j+b+1]);
					const Tad vbottom = 0.25*(vleft+vright+v[i+b-1][j+b-1]+v[i+b][j+b-1]);

					ubar = 0.25*(uleft+uright+utop+ubottom);
					vbar = 0.25*(vleft+vright+vtop+vbottom);

					const Tad Tleft = T[i+b-1][j+b];
					const Tad Tright = T[i+b][j+b];
					const Tad Ttop = 0.25*(Tleft+Tright+T[i+b-1][j+b+1]+T[i+b][j+b+1]);
					const Tad Tbottom = 0.25*(Tleft+Tright+T[i+b-1][j+b-1]+T[i+b][j+b-1]);


						
						double normal_top[2];
						double normal_bottom[2];
						double normal_left[2];
						double normal_right[2];
						// for(int k=0; k<2; k++){
						// 	normal_top[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j+1][k]);
						// 	normal_bottom[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j-1][k]);
						// 	normal_left[k] = 0.5*(mesh->normal_chi[i][j] + mesh->normal_chi[i][j-1]);
						// 	normal_right[k] = 0.5*(mesh->normal_chi[i+1][j] + mesh->normal_chi[i+1][j-1]);
						// }
						double volume = mesh->volume[i-1][j];
						for(int k=0; k<2; k++){
							normal_top[k] = mesh->normal_eta[i-1][j+1][k];
							normal_bottom[k] = mesh->normal_eta[i-1][j][k];
							normal_right[k] = mesh->normal_chi[i][j][k];
							normal_left[k] = 0.5*(mesh->normal_chi[i][j][k] + mesh->normal_chi[i-1][j][k]);
						}
						
						
						dudx = (normal_top[0]*utop - normal_bottom[0]*ubottom + normal_right[0]*uright - normal_left[0]*uleft)/volume;
						dudy = (normal_top[1]*utop - normal_bottom[1]*ubottom + normal_right[1]*uright - normal_left[1]*uleft)/volume;
						
						dvdx = (normal_top[0]*vtop - normal_bottom[0]*vbottom + normal_right[0]*vright - normal_left[0]*vleft)/volume;
						dvdy = (normal_top[1]*vtop - normal_bottom[1]*vbottom + normal_right[1]*vright - normal_left[1]*vleft)/volume;
						
						dTdx = (normal_top[0]*Ttop - normal_bottom[0]*Tbottom + normal_right[0]*Tright - normal_left[0]*Tleft)/volume;
						dTdy = (normal_top[1]*Ttop - normal_bottom[1]*Tbottom + normal_right[1]*Tright - normal_left[1]*Tleft)/volume;

				}
				else{

					/* dvdy = grad_v[i][j][1]; */
					const uint b = 1;

					const Tad uleft = u[i+b-1][j+b];
					const Tad uright = u[i+b][j+b];
					const Tad utop = 0.25*(uleft+uright+u[i+b-1][j+b+1]+u[i+b][j+b+1]);
					const Tad ubottom = 0.25*(uleft+uright+u[i+b-1][j+b-1]+u[i+b][j+b-1]);

					const Tad vleft = v[i+b-1][j+b];
					const Tad vright = v[i+b][j+b];
					const Tad vtop = 0.25*(vleft+vright+v[i+b-1][j+b+1]+v[i+b][j+b+1]);
					const Tad vbottom = 0.25*(vleft+vright+v[i+b-1][j+b-1]+v[i+b][j+b-1]);
					ubar = 0.25*(uleft+uright+utop+ubottom);
					vbar = 0.25*(vleft+vright+vtop+vbottom);

					
					const Tad Tleft = T[i+b-1][j+b];
					const Tad Tright = T[i+b][j+b];
					const Tad Ttop = 0.25*(Tleft+Tright+T[i+b-1][j+b+1]+T[i+b][j+b+1]);
					const Tad Tbottom = 0.25*(Tleft+Tright+T[i+b-1][j+b-1]+T[i+b][j+b-1]);


						
						double normal_top[2];
						double normal_bottom[2];
						double normal_left[2];
						double normal_right[2];
						// for(int k=0; k<2; k++){
						// 	normal_top[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j+1][k]);
						// 	normal_bottom[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j-1][k]);
						// 	normal_left[k] = 0.5*(mesh->normal_chi[i][j] + mesh->normal_chi[i][j-1]);
						// 	normal_right[k] = 0.5*(mesh->normal_chi[i+1][j] + mesh->normal_chi[i+1][j-1]);
						// }
						double volume = 0.5*(mesh->volume[i][j] + mesh->volume[i-1][j]);
						for(int k=0; k<2; k++){
							normal_top[k] = 0.5*(mesh->normal_eta[i][j+1][k] + mesh->normal_eta[i-1][j+1][k]);
							normal_bottom[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i-1][j][k]);
							normal_right[k] = 0.5*(mesh->normal_chi[i][j][k] + mesh->normal_chi[i+1][j][k]);
							normal_left[k] = 0.5*(mesh->normal_chi[i][j][k] + mesh->normal_chi[i-1][j][k]);
						}

						
						dudx = (normal_top[0]*utop - normal_bottom[0]*ubottom + normal_right[0]*uright - normal_left[0]*uleft)/volume;
						dudy = (normal_top[1]*utop - normal_bottom[1]*ubottom + normal_right[1]*uright - normal_left[1]*uleft)/volume;
						
						dvdx = (normal_top[0]*vtop - normal_bottom[0]*vbottom + normal_right[0]*vright - normal_left[0]*vleft)/volume;
						dvdy = (normal_top[1]*vtop - normal_bottom[1]*vbottom + normal_right[1]*vright - normal_left[1]*vleft)/volume;
						
						dTdx = (normal_top[0]*Ttop - normal_bottom[0]*Tbottom + normal_right[0]*Tright - normal_left[0]*Tleft)/volume;
						dTdy = (normal_top[1]*Ttop - normal_bottom[1]*Tbottom + normal_right[1]*Tright - normal_left[1]*Tleft)/volume;

				}

				const Tad tau_xy = mu*(dudy + dvdx);
				const Tad tau_xx = mu*(2.0*dudx - 2.0/3.0*(dudx + dvdy));
				const Tad tau_yy = mu*(2.0*dvdy - 2.0/3.0*(dudx + dvdy));
				const Tad q_x = -k*dTdx;
				const Tad q_y = -k*dTdy;

				flux_xi[i][j][0] = 0.0;
				flux_xi[i][j][1] = tau_xx*nx + tau_xy*ny;
				flux_xi[i][j][2] = tau_xy*nx + tau_yy*ny;
				flux_xi[i][j][3] = (ubar*tau_xx + vbar*tau_xy - q_x)*nx + (ubar*tau_xy + vbar*tau_yy - q_y)*ny;
				//flux_xi[i][j][3] = 0.0;


			}
		}


				
		for(int i=0; i<nic; i++){
			for(int j=0; j<nj; j++){
				Tad dudx, dudy, dvdx, dvdy, ubar, vbar, dTdx, dTdy;
					auto nx = mesh->normal_eta[i][j][0];
					auto ny = mesh->normal_eta[i][j][1];
					
					if(j == 0){
						/* const Tad dudeta = (u[i+1][j+1] - u[i+1][j]);
						const Tad dvdeta = (v[i+1][j+1] - v[i+1][j]);
						const Tad dudxi = (u[i+2][j+1] - u[i+1-1][j+1])/2.0;
						const Tad dvdxi = (v[i+2][j+1] - v[i+1-1][j+1])/2.0;

						*/
						
						//dudx = mesh->xi_x[i][j]*dudxi + mesh->eta_x[i][j]*dudeta;
						//dudy = mesh->xi_y[i][j]*dudxi + mesh->eta_y[i][j]*dudeta;
						
						//dvdx = mesh->xi_x[i][j]*dvdxi + mesh->eta_x[i][j]*dvdeta;
						//dvdy = mesh->xi_y[i][j]*dvdxi + mesh->eta_y[i][j]*dvdeta;
						
						//dudx = 0.0; dudy=0.0; dvdx=0.0; dvdy=0.0;
						
					   /* 	dudx = grad_u[i][j][0]; */
					/* dudy = grad_u[i][j][1]; */
					
					/* dvdx = grad_v[i][j][0]; */
					/* dvdy = grad_v[i][j][1]; */
						const uint b = 1;
						const Tad utop = u[i+b][j+b];
						const Tad ubottom = u[i+b][j+b-1];

						const Tad uleft = 0.25*(utop + ubottom + u[i+b-1][j+b] + u[i+b-1][j+b-1]);
						const Tad uright = 0.25*(utop + ubottom + u[i+b+1][j+b] + u[i+b+1][j+b-1]);

						const Tad vtop = v[i+b][j+b];
						const Tad vbottom = v[i+b][j+b-1];

						const Tad vleft = 0.25*(vtop + vbottom + v[i+b-1][j+b] + v[i+b-1][j+b-1]);
						const Tad vright = 0.25*(vtop + vbottom + v[i+b+1][j+b] + v[i+b+1][j+b-1]);


						const Tad Ttop = T[i+b][j+b];
						const Tad Tbottom = T[i+b][j+b-1];

						const Tad Tleft = 0.25*(Ttop + Tbottom + T[i+b-1][j+b] + T[i+b-1][j+b-1]);
						const Tad Tright = 0.25*(Ttop + Tbottom + T[i+b+1][j+b] + T[i+b+1][j+b-1]);
						ubar = 0.25*(uleft+uright+utop+ubottom);
						vbar = 0.25*(vleft+vright+vtop+vbottom);

						
						
						double normal_top[2];
						double normal_bottom[2];
						double normal_left[2];
						double normal_right[2];
						// for(int k=0; k<2; k++){
						// 	normal_top[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j+1][k]);
						// 	normal_bottom[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j-1][k]);
						// 	normal_left[k] = 0.5*(mesh->normal_chi[i][j] + mesh->normal_chi[i][j-1]);
						// 	normal_right[k] = 0.5*(mesh->normal_chi[i+1][j] + mesh->normal_chi[i+1][j-1]);
						// }
						double volume = mesh->volume[i][j];
						for(int k=0; k<2; k++){
							normal_top[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j+1][k]);
							normal_bottom[k] = mesh->normal_eta[i][j][k];
							normal_left[k] = mesh->normal_chi[i][j][k];
							normal_right[k] = mesh->normal_chi[i+1][j][k];
						}

						
						dudx = (normal_top[0]*utop - normal_bottom[0]*ubottom + normal_right[0]*uright - normal_left[0]*uleft)/volume;
						dudy = (normal_top[1]*utop - normal_bottom[1]*ubottom + normal_right[1]*uright - normal_left[1]*uleft)/volume;
						
						dvdx = (normal_top[0]*vtop - normal_bottom[0]*vbottom + normal_right[0]*vright - normal_left[0]*vleft)/volume;
						dvdy = (normal_top[1]*vtop - normal_bottom[1]*vbottom + normal_right[1]*vright - normal_left[1]*vleft)/volume;

						
						dTdx = (normal_top[0]*Ttop - normal_bottom[0]*Tbottom + normal_right[0]*Tright - normal_left[0]*Tleft)/volume;
						dTdy = (normal_top[1]*Ttop - normal_bottom[1]*Tbottom + normal_right[1]*Tright - normal_left[1]*Tleft)/volume;

						
				}

				else if(j == nj-1){
					
					

					const uint b = 1;
					const Tad utop = u[i+b][j+b];
					const Tad ubottom = u[i+b][j+b-1];

					const Tad uleft = 0.25*(utop + ubottom + u[i+b-1][j+b] + u[i+b-1][j+b-1]);
					const Tad uright = 0.25*(utop + ubottom + u[i+b+1][j+b] + u[i+b+1][j+b-1]);

					const Tad vtop = v[i+b][j+b];
					const Tad vbottom = v[i+b][j+b-1];

					const Tad vleft = 0.25*(vtop + vbottom + v[i+b-1][j+b] + v[i+b-1][j+b-1]);
					const Tad vright = 0.25*(vtop + vbottom + v[i+b+1][j+b] + v[i+b+1][j+b-1]);

					const Tad Ttop = T[i+b][j+b];
					const Tad Tbottom = T[i+b][j+b-1];

					const Tad Tleft = 0.25*(Ttop + Tbottom + T[i+b-1][j+b] + T[i+b-1][j+b-1]);
					const Tad Tright = 0.25*(Ttop + Tbottom + T[i+b+1][j+b] + T[i+b+1][j+b-1]);

					ubar = 0.25*(uleft+uright+utop+ubottom);
					vbar = 0.25*(vleft+vright+vtop+vbottom);

						
					double normal_top[2];
					double normal_bottom[2];
					double normal_left[2];
					double normal_right[2];
					for(int k=0; k<2; k++){
						normal_top[k] = mesh->normal_eta[i][j][k];
						normal_bottom[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j-1][k]);
						normal_left[k] = mesh->normal_chi[i][j-1][k];
						normal_right[k] = mesh->normal_chi[i+1][j-1][k];
					}
					double volume = mesh->volume[i][j-1];
					dudx = (normal_top[0]*utop - normal_bottom[0]*ubottom + normal_right[0]*uright - normal_left[0]*uleft)/volume;
					dudy = (normal_top[1]*utop - normal_bottom[1]*ubottom + normal_right[1]*uright - normal_left[1]*uleft)/volume;
						
					dvdx = (normal_top[0]*vtop - normal_bottom[0]*vbottom + normal_right[0]*vright - normal_left[0]*vleft)/volume;
					dvdy = (normal_top[1]*vtop - normal_bottom[1]*vbottom + normal_right[1]*vright - normal_left[1]*vleft)/volume;
					
					dTdx = (normal_top[0]*Ttop - normal_bottom[0]*Tbottom + normal_right[0]*Tright - normal_left[0]*Tleft)/volume;
					dTdy = (normal_top[1]*Ttop - normal_bottom[1]*Tbottom + normal_right[1]*Tright - normal_left[1]*Tleft)/volume;

				}
				else{


					const uint b = 1;
					const Tad utop = u[i+b][j+b];
					const Tad ubottom = u[i+b][j+b-1];

						const Tad uleft = 0.25*(utop + ubottom + u[i+b-1][j+b] + u[i+b-1][j+b-1]);
						const Tad uright = 0.25*(utop + ubottom + u[i+b+1][j+b] + u[i+b+1][j+b-1]);

						const Tad vtop = v[i+b][j+b];
						const Tad vbottom = v[i+b][j+b-1];

						const Tad vleft = 0.25*(vtop + vbottom + v[i+b-1][j+b] + v[i+b-1][j+b-1]);
						const Tad vright = 0.25*(vtop + vbottom + v[i+b+1][j+b] + v[i+b+1][j+b-1]);

						const Tad Ttop = T[i+b][j+b];
						const Tad Tbottom = T[i+b][j+b-1];

						const Tad Tleft = 0.25*(Ttop + Tbottom + T[i+b-1][j+b] + T[i+b-1][j+b-1]);
						const Tad Tright = 0.25*(Ttop + Tbottom + T[i+b+1][j+b] + T[i+b+1][j+b-1]);
						ubar = 0.25*(uleft+uright+utop+ubottom);
					vbar = 0.25*(vleft+vright+vtop+vbottom);

						
						double normal_top[2];
						double normal_bottom[2];
						double normal_left[2];
						double normal_right[2];
						for(int k=0; k<2; k++){
							normal_top[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j+1][k]);
							normal_bottom[k] = 0.5*(mesh->normal_eta[i][j][k] + mesh->normal_eta[i][j-1][k]);
							normal_left[k] = 0.5*(mesh->normal_chi[i][j][k] + mesh->normal_chi[i][j-1][k]);
							normal_right[k] = 0.5*(mesh->normal_chi[i+1][j][k] + mesh->normal_chi[i+1][j-1][k]);
						}
						double volume = 0.5*(mesh->volume[i][j] + mesh->volume[i][j-1]);
						dudx = (normal_top[0]*utop - normal_bottom[0]*ubottom + normal_right[0]*uright - normal_left[0]*uleft)/volume;
						dudy = (normal_top[1]*utop - normal_bottom[1]*ubottom + normal_right[1]*uright - normal_left[1]*uleft)/volume;
						
						dvdx = (normal_top[0]*vtop - normal_bottom[0]*vbottom + normal_right[0]*vright - normal_left[0]*vleft)/volume;
						dvdy = (normal_top[1]*vtop - normal_bottom[1]*vbottom + normal_right[1]*vright - normal_left[1]*vleft)/volume;

						dTdx = (normal_top[0]*Ttop - normal_bottom[0]*Tbottom + normal_right[0]*Tright - normal_left[0]*Tleft)/volume;
						dTdy = (normal_top[1]*Ttop - normal_bottom[1]*Tbottom + normal_right[1]*Tright - normal_left[1]*Tleft)/volume;

				}

				const Tad tau_xy = mu*(dudy + dvdx);
				const Tad tau_xx = mu*(2.0*dudx - 2.0/3.0*(dudx + dvdy));
				const Tad tau_yy = mu*(2.0*dvdy - 2.0/3.0*(dudx + dvdy));
				const Tad q_x = -k*dTdx;
				const Tad q_y = -k*dTdy;


				flux_eta[i][j][0] = 0.0;
				flux_eta[i][j][1] = tau_xx*nx + tau_xy*ny;
				flux_eta[i][j][2] = tau_xy*nx + tau_yy*ny;
				//spdlog::get("console")->debug("dudy = {}, dvdx = {}, ubar = {} ny = {}", dudy.value(), dvdx.value(), ubar.value(), ny);
				//flux_eta[i][j][3] = (ubar*tau_xx + vbar*tau_xy - q_x)*nx + (ubar*tau_xy + vbar*tau_yy - q_y)*ny;
				flux_eta[i][j][3] = (ubar*tau_xx + vbar*tau_xy - q_x)*nx + (ubar*tau_xy + vbar*tau_yy - q_y)*ny;
				//flux_eta[i][j][3] = 0.0;
			}
		}


		for(uint i=0; i< nic; i++){
			for(uint j=0; j< njc; j++){
				for(uint k=1; k<mesh->solution->nq; k++){
					a_rhs[i][j][k] += (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
					a_rhs[i][j][k] += (flux_xi[i+1][j][k] - flux_xi[i][j][k]);
					if(k==3)
						spdlog::get("console")->debug("rhs = {}, {}, {}",  (flux_eta[i][j+1][k] - flux_eta[i][j][k]).value(), (flux_xi[i+1][j][k] - flux_xi[i][j][k]).value(), mesh->volume[i][j]);
					//spdlog::get("console")->debug("viscous contribution {} {} ", (flux_eta[i][j+1][k].value() - flux_eta[i][j][k].value()), (flux_xi[i+1][j][k].value() - flux_xi[i][j][k].value()));
				}
			}
		}
	};
	void calc_source_residual(Array3D<Tad>& a_rhs){
		for(uint i=0; i< nic; i++){
			for(uint j=0; j< njc; j++){
				const Tad dpdx = -.1;
				a_rhs[i][j][1] += -dpdx*mesh->volume[i][j];
			}
		}
	};
	void calc_primvars(Array3D<Tad>& a_q);
	void calc_boundary();
	EulerEquation(std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
		mesh = val_mesh;
		config = val_config;

		ni = mesh->ni;
		nj = mesh->nj;
		nq = mesh->solution->nq;
		nic = ni - 1;
		njc = nj - 1;
		
		rho = Array2D<Tad>(nic+2, njc+2);
		u = Array2D<Tad>(nic+2, njc+2);
		v = Array2D<Tad>(nic+2, njc+2);
		p = Array2D<Tad>(nic+2, njc+2);
		T = Array2D<Tad>(nic+2, njc+2);
		
		rholft_xi = Array2D<Tad>(ni, njc);
		ulft_xi = Array2D<Tad>(ni, njc);
		vlft_xi = Array2D<Tad>(ni, njc);
		plft_xi = Array2D<Tad>(ni, njc);
		
		rhorht_xi = Array2D<Tad>(ni, njc);
		urht_xi = Array2D<Tad>(ni, njc);
		vrht_xi = Array2D<Tad>(ni, njc);
		prht_xi = Array2D<Tad>(ni, njc);
		
		rholft_eta = Array2D<Tad>(nic, nj);
		ulft_eta = Array2D<Tad>(nic, nj);
		vlft_eta = Array2D<Tad>(nic, nj);
		plft_eta = Array2D<Tad>(nic, nj);
		
		rhorht_eta = Array2D<Tad>(nic, nj);
		urht_eta = Array2D<Tad>(nic, nj);
		vrht_eta = Array2D<Tad>(nic, nj);
		prht_eta = Array2D<Tad>(nic, nj);

		

		grad_u = Array3D<Tad>(nic, njc, 3);
		grad_v = Array3D<Tad>(nic, njc, 3);
		grad_T = Array3D<Tad>(nic, njc, 3);
		

		flux_xi = Array3D<Tad>(ni, njc, 4U);
		flux_eta = Array3D<Tad>(nic, nj, 4U);

		//transport = std::make_unique<TransportEquation<Tx, Tad>>(mesh, config);

		
		if(config->solver->order == 1){
			reconstruction = std::make_unique<FirstOrder<Tx, Tad>>(ni, nj);
		}
		else if(config->solver->order == 2){
			reconstruction = std::make_unique<SecondOrder<Tx, Tad>>(ni, nj);
		}
		else{
			
		}

		if(config->solver->flux == "roe")
			convective_flux = std::make_unique<RoeFlux<Tx, Tad>>();
		else if(config->solver->flux == "ausm")
			convective_flux = std::make_unique<AUSMFlux<Tx, Tad>>();
		else
			spdlog::get("console")->critical("Flux not found.");
		
		for(const auto& bc : config->geometry->boundary){
			auto name = bc->name;
			auto type = bc->type;
			auto face = bc->face;
			auto start = bc->start;
			auto end = bc->end;
			uint facei;
			if(face=="left") facei = left;
			if(face=="right") facei = right;
			if(face=="bottom") facei = bottom;
			if(face=="top") facei = top;
			if(type == "freestream"){
				auto boundarycondition = new BoundaryConditionFreestream<Tx, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else if(type=="outflow"){
				auto boundarycondition = new BoundaryConditionOutflow<Tx, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else if(type=="periodic"){
				auto boundarycondition = new BoundaryConditionPeriodic<Tx, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else if(type == "slipwall"){
				auto boundarycondition = new BoundaryConditionInviscidWall<Tx, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else if(type == "wall"){
				auto boundarycondition = new BoundaryConditionViscousWall<Tx, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else if(type == "wake"){
				auto boundarycondition = new BoundaryConditionWake<Tx, Tad>(name, mesh, config, facei, start, end);
				boundaryconditions.push_back(boundarycondition);
			}
			else{
				spdlog::get("console")->info("Wrong type of BC.");
			}
		}
		
	};

	~EulerEquation(){
		for(auto&& bc : boundaryconditions)
			delete bc;
	};
};
template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_convective_residual(Array3D<Tad>& a_rhs){
	convective_flux->evaluate(mesh->normal_chi,
							  rholft_xi, ulft_xi, vlft_xi, plft_xi,
							  rhorht_xi, urht_xi, vrht_xi, prht_xi,
							  flux_xi);
	convective_flux->evaluate(mesh->normal_eta,
							  rholft_eta, ulft_eta, vlft_eta, plft_eta,
							  rhorht_eta, urht_eta, vrht_eta, prht_eta,
							  flux_eta);
#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq; k++){
				a_rhs[i][j][k] -= (flux_eta[i][j+1][k] - flux_eta[i][j][k]);
				a_rhs[i][j][k] -= (flux_xi[i+1][j][k] - flux_xi[i][j][k]);
			}
		}
	}
}

template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_primvars(Array3D<Tad>& a_q){
	primvars<Tad>(a_q, rho, u, v, p, T, 1U, 1U);
}


template <class Tx, class Tad>
void EulerEquation<Tx, Tad>::calc_boundary(){
	for(auto&& bc : boundaryconditions)
		bc->apply(rho, u, v, p, T);
}

template <class Tx, class Tad>
	void EulerEquation<Tx, Tad>::calc_residual(Array3D<Tad>& a_q, Array3D<Tad>& a_rhs){
	a_rhs.fill(0.0);

	calc_primvars(a_q);
	calc_boundary();
	
	reconstruction->evaluate_xi(rho, rholft_xi, rhorht_xi);
	reconstruction->evaluate_xi(u, ulft_xi, urht_xi);
	reconstruction->evaluate_xi(v, vlft_xi, vrht_xi);
	reconstruction->evaluate_xi(p, plft_xi, prht_xi);

	reconstruction->evaluate_eta(rho, rholft_eta, rhorht_eta);
	reconstruction->evaluate_eta(u, ulft_eta, urht_eta);
	reconstruction->evaluate_eta(v, vlft_eta, vrht_eta);
	reconstruction->evaluate_eta(p, plft_eta, prht_eta);
	
	calc_convective_residual(a_rhs);
	calc_viscous_residual(a_rhs);
	calc_source_residual(a_rhs);
	//transport->calc_residual(a_q, a_rhs, u, v,
	//						 ulft_xi, urht_xi,
	//						 ulft_eta, urht_eta,
	//						 vlft_xi, vrht_xi,
	//						 vlft_eta, vrht_eta);
#pragma omp parallel for
	for(uint i=0; i< nic; i++){
		for(uint j=0; j< njc; j++){
			for(uint k=0; k<mesh->solution->nq+mesh->solution->ntrans; k++){
				a_rhs[i][j][k] /= mesh->volume[i][j];
			}
		}
	}

};

#endif
