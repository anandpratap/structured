#ifndef _IO_H
#define _IO_H
#include "common.h"
#include "mesh.h"
#include "config.h"

template<class Tx>
class IOManager{
public:
	std::shared_ptr<Config<Tx>> config;
	std::shared_ptr<Mesh<Tx>> mesh;
	std::string label;
	IOManager(std::shared_ptr<Mesh<Tx>> val_mesh, std::shared_ptr<Config<Tx>> val_config){
		config = val_config;
		mesh = val_mesh;

		// required for npz
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ntrans = mesh->solution->ntrans;

		if(mesh->label == ""){
			label = config->io->label;
		}
		else{
			label = config->io->label + "_" + mesh->label;
		}
	};

	~IOManager(){
		
	}

	void write(const uint iteration){
		write_tecplot();
		write_npz();
		write_restart();
		write_surface();
	}
	
	void write_tecplot(){
		std::string filename = label + ".tec";
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		auto q = mesh->solution->q;
		auto xc = mesh->xc;
		auto yc = mesh->yc;

		const auto nq = mesh->solution->nq;
		const auto ntrans = mesh->solution->ntrans;

		std::ofstream outfile;
		outfile.open(filename);
		char buffer [500];
		outfile<<"title = \"Solution\""<<"\n";
		outfile<<"variables = \"x\" \"y\"";
		for(uint tn=0; tn<nq; tn++){
			outfile << " \""<<mesh->solution->q_name[tn]<<"\"";
		}
		for(uint tn=0; tn<ntrans; tn++){
			outfile << " \"psi_0\"";
		}
		for(uint tn=0; tn<mesh->solution->naux; tn++){
			outfile << " \""<<mesh->solution->q_aux_name[tn]<<"\"";
		}

		outfile << "\n";

		outfile<<"zone i="<<nic<<", j="<<njc<<", f=point\n";
		
		for(uint j=0; j<njc; j++){
			for(uint i=0; i<nic; i++){
				outfile << xc[i][j] << " ";
				outfile << yc[i][j] << " ";

				for(uint k=0; k<nq; k++){
					outfile << q[i][j][k] << " ";
				}

				for(uint tn=0; tn<ntrans; tn++){
					outfile << q[i][j][4+tn] << " ";
				}

				for(uint tn=0; tn<mesh->solution->naux; tn++){
					outfile << mesh->solution->q_aux[i][j][tn] << " ";
				}
		
				outfile << "\n";
			}
		}
		
		outfile.close();
		spdlog::get("console")->info("Wrote tecplot file {}.", filename);
	}
	
	
	void write_npz(){
		std::string filename = label + ".npz";
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		auto q = mesh->solution->q;
		auto xc = mesh->xc;
		auto yc = mesh->yc;
		const unsigned int shape[] = {nic, njc};
		const unsigned int shapeq[] = {nic, njc, 4};
		auto rho = mesh->solution->rho;
		auto u = mesh->solution->u;
		auto v = mesh->solution->v;
		auto p = mesh->solution->p;
		auto T = mesh->solution->T;
		
		cnpy::npz_save(filename,"xc",xc.data(),shape,2,"w");
		cnpy::npz_save(filename,"yc",yc.data(),shape,2,"a");
		cnpy::npz_save(filename,"q",q.data(),shapeq,3,"a");
		cnpy::npz_save(filename,"rho",rho.data(),shape,2,"a");
		cnpy::npz_save(filename,"u",u.data(),shape,2,"a");
		cnpy::npz_save(filename,"v",v.data(),shape,2,"a");
		cnpy::npz_save(filename,"p",p.data(),shape,2,"a");
		cnpy::npz_save(filename,"T",T.data(),shape,2,"a");
		
		spdlog::get("console")->info("Wrote numpy npz file {}.", filename);
	}


	void write_restart(){
		std::string filename = label + ".out";
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ntrans = mesh->solution->ntrans;
		auto q = mesh->solution->q;
		auto xc = mesh->xc;
		auto yc = mesh->yc;
		std::ofstream outfile(filename,std::ofstream::binary);
		for(int i=0; i<nic; i++){
			for(int j=0; j<njc; j++){
				outfile.write(reinterpret_cast<const char*>(&q[i][j][0]), sizeof(Tx)*(nq+ntrans));
			}
		}
		outfile.close();
		spdlog::get("console")->info("Wrote restart file {}.", filename);
	}

	void read_restart(){
		std::string filename = label + ".out";
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		const auto ntrans = mesh->solution->ntrans;
				
		auto q = mesh->solution->q;
		auto xc = mesh->xc;
		auto yc = mesh->yc;
		std::ifstream infile(filename,std::ofstream::binary);

		infile.seekg (0,infile.end);
		long size = infile.tellg();
		long size_expected = nic*njc*(nq+ntrans)*sizeof(Tx);
		infile.seekg (0);
		assert(size == size_expected);
		
		for(int i=0; i<nic; i++){
			for(int j=0; j<njc; j++){
				infile.read(reinterpret_cast<char*>(&q[i][j][0]), sizeof(Tx)*(nq+ntrans));
			}
		}
		infile.close();
		
	}
	void write_surface(){
		std::string filename = label + ".surface";
		const auto nic = mesh->nic;
		const auto njc = mesh->njc;
		const auto nq = mesh->solution->nq;
		auto q = mesh->solution->q;
		auto xc = mesh->xc;
		auto yc = mesh->yc;
		Tx p_inf = 1/1.4;
		Tx rho_inf = 1.0;
		Tx u_inf = 0.5;
		int j1 = mesh->j1-1;
		Tx x, cp;
		std::ofstream outfile;
		outfile.open(filename);
		auto rho = mesh->solution->rho;
		auto u = mesh->solution->u;
		auto v = mesh->solution->v;
		auto p = mesh->solution->p;
		auto T = mesh->solution->T;
		for(uint i=j1; i<j1+mesh->nb; i++){
			x = xc[i][0];
			cp = (p[i][0] - p_inf)/(0.5*rho_inf*u_inf*u_inf);
			outfile<<x<<" "<<cp<<std::endl;
		}
		outfile.close();
		spdlog::get("console")->info("Wrote surface file {}.", filename);
	}
};

#endif
