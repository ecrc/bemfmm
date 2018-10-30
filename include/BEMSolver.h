#ifndef BEMSOLVER
#define BEMSOLVER
#include "utils.h"
#include "fbind.h"
#include "gmres.h"
#include "fmm/acoustics_wrapper.h"
#include <fstream>
#include "global_data.h"

namespace bemfmm {
	enum field_mode {
		FAR = 0,
		NEAR_NEAR_SCHEME = 1,
		NEAR_FAR_SCHEME = 2
	};
	class BEMSolver {
	private:
		integral_data const& int_data;
		common_data & comm_data;     
		int** const& nsupan;
		double** const& sunod;
		int const& ntriangle;
		int const& nnpt;
		int const& nipp;
		d_complex_t const& cjvk;
		double const& nearpd;
		double vkk[3];
		double vhh[3];
		int const& nhdgqp; 
		int const& nlqp; 
		int const& mpirank;
		int const& mpisize;			
		int_data_type const& alphao;
		int_data_type const& betao;
		int_data_type const& gammao;
		int_data_type const& alphas;
		int_data_type const& betas;
		int_data_type const& gammas;
		int_data_type const& alphas1;
		int_data_type const& betas1;
		int_data_type const& gammas1;			
		d_vector const& wwo;
		d_vector const& xxo;
		d_vector const& wws;
		d_vector const& xxs; 
		d_vector const& ws;
    d_vector const& ws1; 
		d_vector xb;
		d_vector yb;
		d_vector zb;
		d_complex_t_vec zzsparse;
		d_complex_t_vec self_contributions;
		d_complex_t_vec rj;
#if USE_PART
		void setRHS(d_complex_t_vec& crhs, int const& partition_size, d_vector const& xt, d_vector const& yt, d_vector const& zt, int16_vec const& pointlocs) {
			d_complex_t I(0.,1.);
			double vos1[2] = {1.0,0.0};
			double vos2[2] = {0.0,1.0};
			double vos3[2] = {0.0,0.0};
			double vippo[3] = {0.0,0.0,0.0};
			double vcsio[2];
			crhs.resize(partition_size,0.0);
			double dot_product;
			for (int i = 0; i < partition_size; ++i) {
				int j = pointlocs[i];
				vcsio[0] = alphao[j]*vos1[0] + betao[j]*vos2[0] + gammao[j]*vos3[0];
				vcsio[1] = alphao[j]*vos1[1] + betao[j]*vos2[1] + gammao[j]*vos3[1];
				utils::getR(vcsio, xt, yt, zt, i/nnpt, nnpt, vippo);
				utils::dot_product(vkk, vippo, dot_product);
				crhs[i]=exp(I*cjvk*dot_product);
			}
		}
#else  	
		void setRHS(d_complex_t_vec& crhs, int const& partition_size, int32_vec const& patches, int16_vec const& pointlocs) {
			d_complex_t I(0.,1.);
			double vos1[2] = {1.0,0.0};
			double vos2[2] = {0.0,1.0};
			double vos3[2] = {0.0,0.0};
			double vippo[3] = {0.0,0.0,0.0};
			double vcsio[2];
			int* nodes_sn;
			crhs.resize(partition_size,0.0);
			double dot_product;
			for (int i = 0; i < partition_size; ++i) {
				nodes_sn = nsupan[patches[i]];
				int j = pointlocs[i];
				vcsio[0] = alphao[j]*vos1[0] + betao[j]*vos2[0] + gammao[j]*vos3[0];
				vcsio[1] = alphao[j]*vos1[1] + betao[j]*vos2[1] + gammao[j]*vos3[1];
				utils::getR(vcsio, nodes_sn, comm_data.sunod,  nnpt, vippo);
				utils::dot_product(vkk, vippo, dot_product);
				crhs[i]=exp(I*cjvk*dot_product);
			}
		}
#endif

		void setNearQuadPoints(d_vector_2d& gaussPoints, d_vector_2d& ipolator_near, double* ipolymatrix, int ntriangle, int const& nipp) {      
			int nhdgqp = comm_data.nhdgqp;
			ipolator_near.resize(nipp,d_vector(nhdgqp,0));
			d_vector_2d vcsik(nhdgqp,d_vector(2,0.0));
			gaussPoints.resize(ntriangle*nhdgqp,std::vector<double>(3));
			utils::initBasisVector(nhdgqp, nipp, int_data.alphas1, int_data.betas1, int_data.gammas1, ipolymatrix, ipolator_near, vcsik);
#pragma omp parallel for 
			for (int i = 0; i < ntriangle; ++i) {
#if USE_PART
				double* _sunod = comm_data.triangles + (i*nipp*3);
#endif        
				for (int m = 0; m < nhdgqp; ++m) {
#if USE_PART
					utils::getR(vcsik[m], _sunod,  nipp, gaussPoints[i*nhdgqp+m]);                                   
#else 
					utils::getR(vcsik[m], nsupan[i], comm_data.sunod,  nipp, gaussPoints[i*nhdgqp+m]);                                   
#endif
				}                  
			}
		}

	public:
		BEMSolver(integral_data const& _int_data, common_data & _comm_data) 
		: int_data(_int_data), comm_data(_comm_data), nsupan(comm_data.nsupan), 
		sunod(comm_data.sunod), ntriangle(comm_data.ntriangle), nnpt(comm_data.nnpt), 
		nipp(comm_data.nipp), cjvk(comm_data.cjvk), nearpd(comm_data.nearpd), 
		nhdgqp(comm_data.nhdgqp), nlqp(comm_data.nlqp), mpirank(comm_data.mpirank), 
		mpisize(comm_data.mpisize), alphao(int_data.alphao), betao(int_data.betao), gammao(int_data.gammao), 
		alphas(int_data.alphas), betas(int_data.betas), gammas(int_data.gammas), 
		alphas1(int_data.alphas1), betas1(int_data.betas1), gammas1(int_data.gammas1),
		wwo(int_data.wwo), xxo(int_data.xxo), wws(int_data.wws), xxs(int_data.xxs), 
		ws(int_data.ws), ws1(int_data.ws1)
		{
#if USE_PART
			size_t const n = nipp * init_part_size;  
			xb.resize(n); yb.resize(n); zb.resize(n);
			utils::fetchCoords(ntriangle, nipp, nnpt, 
				alphas, betas, gammas, comm_data, xb, yb, zb);   
#else
			size_t const n = nipp * ntriangle;  
			xb.resize(n); yb.resize(n); zb.resize(n);
			utils::fetchCoords(ntriangle, nipp, nnpt, nsupan, 
				alphas, betas, gammas, sunod, xb, yb, zb); 
#endif
		  for (int i = 0; i < 3; ++i) {
        vkk[i] = comm_data.vkk[i];
        vhh[i] = comm_data.vhh[i];
      }
		}

		void buildDenseMatrix() {
			zzsparse.resize((nipp*ntriangle) * (nipp*ntriangle));          
			self_contributions.resize(nipp * nipp * ntriangle);
			d_ptr ruv1, ruv2;
			d_complex_t I(0.,1.);
			double veuv1[2]={0.0, 0.0},veuv2[2]={0.0, 0.0},jacob_star;
			double vcsik[2]={0.0, 0.0};      
			int* nodeo_sn;
			int* nodes_sn;
			double vippo[3]={0.0, 0.0, 0.0},vipps[3]={0.0, 0.0, 0.0},vgaussk[3]={0.0, 0.0, 0.0};
			double vuo[3]={0.0, 0.0, 0.0},vvo[3]={0.0, 0.0, 0.0};
			double vus[3]={0.0, 0.0, 0.0};
			double vvs[3]={0.0, 0.0, 0.0};
			double vrrr[3]={0.0, 0.0, 0.0}, rrr;
			int ielement00;
			d_complex_t z00;
			double polyvector[nipp];
			double ipolator; 
			double ppp,qqq,aa1[3]={0.0, 0.0, 0.0},aa2[3]={0.0, 0.0, 0.0},bb1[3]={0.0, 0.0, 0.0},vmmm[3]={0.0, 0.0, 0.0},valf[3]={0.0, 0.0, 0.0};
			double vrr0[3]={0.0, 0.0, 0.0},rr0;
			double polymatrix[nipp*nipp];
			double ipolymatrix[nipp*nipp]; 
			d_vector_2d ipolator_near(nipp,d_vector(nhdgqp,0));
			d_vector_2d vcsik_(nhdgqp,d_vector(2,0.0));
			d_vector_2d vcsic(nhdgqp,d_vector(2,0.0));
			d_vector_2d vcsis(nipp,d_vector(2,0.0));
			d_vector_2d vcsio(nipp,d_vector(2,0.0));
			for (int k = 0; k < nipp; ++k) {
				vcsio[k][0] = alphao[k] * utils::vos1[0] + betao[k] * utils::vos2[0] + gammao[k] * utils::vos3[0];
				vcsio[k][1] = alphao[k] * utils::vos1[1] + betao[k] * utils::vos2[1] + gammao[k] * utils::vos3[1];
			}
			for (int k = 0; k < nipp; ++k) {
				vcsis[k][0] = alphas[k] * utils::vos1[0] + betas[k] * utils::vos2[0] + gammas[k] * utils::vos3[0];
				vcsis[k][1] = alphas[k] * utils::vos1[1] + betas[k] * utils::vos2[1] + gammas[k] * utils::vos3[1];
			}
			utils::initPolyMatrix(polymatrix, nipp, alphas, betas, gammas, ipolymatrix);
			utils::initBasisVector(nhdgqp, nipp, alphas1, betas1, gammas1, ipolymatrix, ipolator_near, vcsik_);    
			for (int i = 0; i < ntriangle; ++i) {
				nodeo_sn = nsupan[i];
#pragma omp parallel for default(none) shared(vus, ipolymatrix, i, I, ipolator_near, utils::vos1, utils::vos2, utils::vos3, vmmm) firstprivate(nodes_sn, vippo, nodeo_sn, vcsio, vipps, vcsis, vvs, vrrr, rrr, ielement00, z00, vcsik_, vgaussk, aa1, aa2, bb1, vvo, vuo, ruv1, ruv2, veuv1, veuv2, jacob_star, ppp, qqq, vcsik, vrr0, rr0, ipolator, polyvector, valf)
				for (int j = 0; j < ntriangle; ++j) {
					nodes_sn = nsupan[j];
					if(i!=j) {        
						for (int k = 0; k < nipp; ++k) {
#if !USE_PART            
							utils::getR(vcsio[k], nodeo_sn, sunod,  nnpt, vippo);
#endif

							for (int l = 0; l < nipp; ++l) {
#if !USE_PART
								utils::getR(vcsis[l], nodes_sn, sunod,  nnpt, vipps);
#endif
								utils::getUV_vector(vcsis[l], nodes_sn, sunod, nnpt, vus, vvs);
								vrrr[0]=vippo[0]-vipps[0];
								vrrr[1]=vippo[1]-vipps[1];
								vrrr[2]=vippo[2]-vipps[2];
								utils::rootDotProduct(vrrr,vrrr,rrr);
								ielement00 = (i * nipp + k) * ntriangle * nipp + j * nipp + l;
								z00 = 0.0;
								if(rrr > nearpd) {
                    z00 = ws[l] * 0.5 * exp(I*cjvk * rrr)/(4.0 * M_PI * rrr); //matrix fill for far intraction
                  } else {
                  	for (int m = 0; m < nhdgqp; ++m) {
#if !USE_PART
                  		utils::getR(vcsik_[m], nodes_sn, sunod,  nnpt, vgaussk);
#endif
                  		vrrr[0] = vippo[0] - vgaussk[0];
                  		vrrr[1] = vippo[1] - vgaussk[1];
                  		vrrr[2] = vippo[2] - vgaussk[2];
                  		utils::rootDotProduct(vrrr,vrrr,rrr);
                  		z00 += ws1[m] * 0.5 * std::exp(I*cjvk * rrr)/(4.0 * M_PI * rrr) * ipolator_near[l][m];
                  	}                    
                  }
                  zzsparse[ielement00] += z00;
                }
              }
            } else { // self-interaction
              //continue;
            	for (int d = 0; d < 3; ++d) {
            		aa1[d]=(sunod[d][nodeo_sn[0]]+sunod[d][nodeo_sn[2]]-2.0*sunod[d][nodeo_sn[5]]);
            		aa2[d]=(sunod[d][nodeo_sn[2]]+sunod[d][nodeo_sn[3]]-sunod[d][nodeo_sn[4]]-sunod[d][nodeo_sn[5]]);
            		bb1[d]=(sunod[d][nodeo_sn[1]]+sunod[d][nodeo_sn[2]]-2.0*sunod[d][nodeo_sn[4]]);
            	}
            	for (int k = 0; k < nipp; ++k) {
#if !USE_PART
            		utils::getR(vcsio[k], nodeo_sn, sunod,  nnpt, vippo);
#endif
            		utils::getUV_vector(vcsio[k], nodeo_sn, sunod, nnpt, vuo, vvo);
            		for (int l = 0; l < nipp; ++l) {
            			ielement00 = (i * nipp + k) * ntriangle * nipp + j * nipp + l;
            			z00 = 0.0;
            			for (int d = 0; d < 3; ++d){
            				switch(d) {
            					case 0:
            					ruv1 = utils::vos2;
            					ruv2 = utils::vos3;
            					break;
            					case 1:
            					ruv1 = utils::vos3;
            					ruv2 = utils::vos1;
            					break;
            					case 2:
            					ruv1 = utils::vos1;
            					ruv2 = utils::vos2;
            					break;
            				}
            				veuv1[0] = vcsio[k][0]- ruv2[0];
            				veuv1[1] = vcsio[k][1]- ruv2[1];
            				veuv2[0] = ruv1[0] - ruv2[0];
            				veuv2[1] = ruv1[1] - ruv2[1];
            				jacob_star = std::abs(veuv1[0]*veuv2[1]-veuv1[1]*veuv2[0]);
            				for (int kk = 0; kk < nlqp; ++kk) {
            					ppp = ruv1[0]-vcsio[k][0]+(ruv2[0]-ruv1[0])*xxo[kk];
            					qqq = ruv1[1]-vcsio[k][1]+(ruv2[1]-ruv1[1])*xxo[kk];                         
            					for (int d = 0; d < 3; ++d) {
            						vmmm[d] = 2.0 * aa1[d] * ppp * ppp + 2.0 * bb1[d] * qqq * qqq + 4.0 * aa2[d] * ppp * qqq;
            						valf[d] = vuo[d] * ppp + vvo[d] * qqq;
            					}
            					for (int mm = 0; mm < nlqp; ++mm) {
            						vcsik[0]=(ruv1[0]-vcsio[k][0])*xxs[mm]+(ruv2[0]-ruv1[0])*xxo[kk]*xxs[mm]+vcsio[k][0];
            						vcsik[1]=(ruv1[1]-vcsio[k][1])*xxs[mm]+(ruv2[1]-ruv1[1])*xxo[kk]*xxs[mm]+vcsio[k][1];
            						for (int d = 0; d < 3; ++d) {
            							vrr0[d]=-xxs[mm] * vmmm[d] - valf[d];
            							vrrr[d]= xxs[mm] * vrr0[d];
            						}                      
            						utils::rootDotProduct(vrrr,vrrr,rrr);
            						utils::rootDotProduct(vrr0,vrr0,rr0);                      
            						ipolator = 0.0;
            						if(nipp == 6) {
            							polyvector[0] = 1.0;
            							polyvector[1] = vcsik[0];
            							polyvector[2] = vcsik[1];
            							polyvector[3] = vcsik[0] * vcsik[1];
            							polyvector[4] = vcsik[0] * vcsik[0];
            							polyvector[5] = vcsik[1] * vcsik[1];
            							for (int kkk = 0; kkk < 6; ++kkk) {
            								ipolator+=ipolymatrix[kkk*nipp+l]*polyvector[kkk];
            							}
            							z00 += wwo[kk]*wws[mm]*exp(I*cjvk*rrr)/(4.0*M_PI*rr0)*ipolator*jacob_star;
            						}
            					}
            				}                  
            			}
            			zzsparse[ielement00]+=z00; 
            			self_contributions[(i*nipp*nipp) + (k * nipp) + l] += z00;
            		}
            	}
            }
          }
        }   
      }

      void calculateScatteredField(int mode, int partition_size, bool write_output, int32_vec const& patches, int16_vec const& pointlocs, std::string file_name) {     
      	int ntheta, nphi, ngprcs=0;
      	double the_sta_degree, the_end_degree, phi_sta_degree, phi_end_degree, rnfobp = 0;
      	int_data_type ws, local_alphas, local_betas, local_gammas;      
      	double vcsis[2], vgaussk[3];
      	double dthe,dphi,phi_degree,the_degree;
      	double vekk[3]={0.0, 0.0, 0.0}, vobp[3]={0.0, 0.0, 0.0}, vipps[3]={0.0, 0.0, 0.0},vrrr[3]={0.0, 0.0, 0.0};
      	d_complex_t scattered_f;
      	double polymatrix[nipp*nipp];
      	double ipolymatrix[nipp*nipp];    
      	d_vector_2d ipolator, vcsik;
      	uint16_t const round_off_decimals = 100;
      	double root_dot_pro;
      	double pro;
      	d_complex_t I(0.,1.);      
      	switch (mode) {
      		case NEAR_FAR_SCHEME:
      		ntheta = comm_data.nthe_nf;
      		nphi   = comm_data.nphi_nf;
      		the_sta_degree = comm_data.the_sta_degree_nf;
      		the_end_degree = comm_data.the_end_degree_nf;
      		phi_sta_degree = comm_data.phi_sta_degree_nf;
      		phi_end_degree = comm_data.phi_end_degree_nf;
      		rnfobp = comm_data.rnfobp;
      		ws    = int_data.ws;
      		break;
      		case FAR:        
      		ntheta = comm_data.nthe_rcs;
      		nphi   = comm_data.nphi_rcs;
      		ngprcs = comm_data.ngprcs;
      		ipolator = d_vector_2d(nipp,d_vector(ngprcs,0));
      		vcsik = d_vector_2d(ngprcs,d_vector(2,0));
      		the_sta_degree = comm_data.the_sta_degree;
      		the_end_degree = comm_data.the_end_degree;
      		phi_sta_degree = comm_data.phi_sta_degree;
      		phi_end_degree = comm_data.phi_end_degree;
      		local_alphas = int_data.alphas2;
      		local_betas  = int_data.betas2;
      		local_gammas = int_data.gammas2;
      		ws     = int_data.ws2;
      		break;
      		case NEAR_NEAR_SCHEME:
      		ntheta = comm_data.nthe_nf;
      		nphi   = comm_data.nphi_nf;
      		ngprcs = comm_data.nhdgqp;
      		rnfobp = comm_data.rnfobp;
      		ipolator = d_vector_2d(nipp,d_vector(ngprcs,0));
      		vcsik    = d_vector_2d(ngprcs,d_vector(2,0));
      		the_sta_degree = comm_data.the_sta_degree_nf;
      		the_end_degree = comm_data.the_end_degree_nf;
      		phi_sta_degree = comm_data.phi_sta_degree_nf;
      		phi_end_degree = comm_data.phi_end_degree_nf;
      		local_alphas = int_data.alphas1;
      		local_betas  = int_data.betas1;
      		local_gammas = int_data.gammas1;
      		ws     = int_data.ws1;
      		break;
      		default:
      		std::cout<<"Error: unknow scattered field mode" << std::endl;
      		return;
      		break;
      	}      
      	std::ofstream output_file;
      	if(write_output) output_file.open(utils::createProcessFile(file_name.c_str(),mpirank).c_str());
      	if(mode != NEAR_FAR_SCHEME) {
      		utils::initPolyMatrix(polymatrix, nipp, alphas,betas, gammas, ipolymatrix);
      		utils::initBasisVector(ngprcs, nipp, local_alphas, local_betas, local_gammas, ipolymatrix, ipolator, vcsik);
      	} 
      	dthe = (ntheta == 1)? 0.0 : ROUND_OFF((the_end_degree - the_sta_degree) / ntheta, round_off_decimals);
      	dphi = (  nphi == 1)? 0.0 : ROUND_OFF((phi_end_degree - phi_sta_degree) / nphi, round_off_decimals);
      	phi_degree = phi_sta_degree;
      	for (int i = 0; i < nphi; ++i, phi_degree+=dphi) {
      		the_degree = the_sta_degree;
      		for (int j = 0; j < ntheta; ++j, the_degree+=dthe) {
	          // Compute r,phi,and theta directed unit vectors. 
      			vekk[0] = sin(the_degree * M_PI / 180.0) * cos( phi_degree * M_PI / 180.0);
      			vekk[1] = sin(the_degree * M_PI / 180.0) * sin( phi_degree * M_PI / 180.0);															 
      			vekk[2] = cos(the_degree * M_PI / 180.0);								
      			if(mode == NEAR_FAR_SCHEME || mode == NEAR_NEAR_SCHEME) {
      				vobp[0] = rnfobp * vekk[0];
      				vobp[1] = rnfobp * vekk[1];
      				vobp[2] = rnfobp * vekk[2];
      			}
      			scattered_f = 0.0;
      			double real_scattered = 0.0;
      			double img_scattered  = 0.0;           
	#pragma omp parallel for private(vcsis,vrrr,root_dot_pro,pro,vgaussk,vipps), reduction(+:real_scattered), reduction(+:img_scattered)
      			for (int k = 0; k < partition_size; ++k) {
      				d_complex_t scat = 0.0;
	#if !USE_PART           
      				int* nodes_sn = nsupan[patches[k]];
	#endif            
      				int l = pointlocs[k];            
      				if(mode == NEAR_FAR_SCHEME) {
      					vcsis[0] = alphas[l] * utils::vos1[0] + betas[l] * utils::vos2[0] + gammas[l] * utils::vos3[0];
      					vcsis[1] = alphas[l] * utils::vos1[1] + betas[l] * utils::vos2[1] + gammas[l] * utils::vos3[1];
	#if USE_PART
      					utils::getR(vcsis, comm_data.xt, comm_data.yt, comm_data.zt, k/nnpt,  nnpt, vipps);
	#else
      					utils::getR(vcsis, nodes_sn, comm_data.sunod,  nnpt, vipps);
	#endif
      					vrrr[0] = vobp[0] - vipps[0];
      					vrrr[1] = vobp[1] - vipps[1];
      					vrrr[2] = vobp[2] - vipps[2];                
      					utils::rootDotProduct(vrrr,vrrr,root_dot_pro);              
      					scat=-ws[l] * exp(I * comm_data.cjvk * root_dot_pro) * rj[k] * 0.5 / (4 * M_PI * root_dot_pro);
      				} else if (mode == FAR || mode == NEAR_NEAR_SCHEME) {
      					for (int m = 0; m < ngprcs; ++m) {
	#if USE_PART
      						utils::getR(vcsik[m], comm_data.xt, comm_data.yt,comm_data.zt, k/nnpt,  nnpt, vgaussk);   
	#else
      						utils::getR(vcsik[m], nodes_sn, comm_data.sunod,  nnpt, vgaussk);                  
	#endif 

      						if(mode == NEAR_NEAR_SCHEME) {
      							vrrr[0] = vobp[0] - vgaussk[0];
      							vrrr[1] = vobp[1] - vgaussk[1];
      							vrrr[2] = vobp[2] - vgaussk[2];
      							utils::rootDotProduct(vrrr,vrrr,pro);
      						} else {
      							utils::dot_product(vekk,vgaussk,pro);                    
      						}
      						if(mode == NEAR_NEAR_SCHEME) {          
      							scat +=-ws[m] * exp(I * comm_data.cjvk * pro) * rj[k] * 0.5 / (4.0 * M_PI * pro) * ipolator[l][m];
      						} else {
      							scat +=-ws[m] * exp(I * -comm_data.cjvk * pro) * rj[k] * 0.5 / (4.0 * M_PI ) * ipolator[l][m];
      						}
      					}
      				}
      				real_scattered += std::real(scat);
      				img_scattered  += std::imag(scat);              
      			}
      			scattered_f = d_complex_t(real_scattered, img_scattered);
      			if(write_output)  output_file << std::setprecision(alogger::precis) << the_degree << " " 
      				<< phi_degree << " " << std::real(scattered_f)     << " " 
      			<< std::imag(scattered_f) << " " << std::abs(scattered_f) << std::endl;          
      		}
      	}
      }

#if USE_PART
      d_complex_t_vec computeSourceField(int& out_size, int32_vec& patches, int16_vec& pointlocs) {
      	int32_t size = comm_data.init_part_size*nipp;
      	d_complex_t_vec crhs;
      	out_size = size;      
      	if(size > 0){
      		patches.resize(size);
      		pointlocs.resize(size);
      	}

      	int offset = (ntriangle/comm_data.mpisize)*comm_data.mpirank;
      	for (int i = 0; i < comm_data.init_part_size; ++i) {        
      		for (int j = 0; j < nipp; ++j) {
      			patches[i*nipp + j] = offset + i;
      			pointlocs[i*nipp + j]  = j;
      		}
      	}
#if USE_FMM
      	double eps2 = 0.0;
      	d_vector const& ws = int_data.ws;
      	d_vector_2d vcsio(nipp,d_vector(2,0.0));
      	for (int k = 0; k < nipp; ++k) {
      		vcsio[k][0] = alphao[k] * utils::vos1[0] + betao[k] * utils::vos2[0] + gammao[k] * utils::vos3[0];
      		vcsio[k][1] = alphao[k] * utils::vos1[1] + betao[k] * utils::vos2[1] + gammao[k] * utils::vos3[1];
      	}
      	d_vector_2d gaussPoints;
      	d_vector_2d ipolator_near;
      	double ipolymatrix[nipp*nipp]; 
      	double polymatrix[nipp*nipp];
      	utils::initPolyMatrix(polymatrix, nipp, int_data.alphas, int_data.betas, int_data.gammas, ipolymatrix);      
      	alogger::startTimer("Set Near Quad Points");
      	setNearQuadPoints(gaussPoints, ipolator_near, ipolymatrix, comm_data.init_part_size, nipp);
      	alogger::stopTimer("Set Near Quad Points");
      	num_threads(comm_data.fmmAttributes.nthreads);
      	alogger::startTimer("FMM Initialization");
      	FMM_Init2(eps2, comm_data.fmmAttributes, size, 
      		DOUBLE_ADDRESS_OFF(xb), DOUBLE_ADDRESS_OFF(yb), DOUBLE_ADDRESS_OFF(zb),  
      		DOUBLE_ADDRESS_OFF(comm_data.xt), DOUBLE_ADDRESS_OFF(comm_data.yt), DOUBLE_ADDRESS_OFF(comm_data.zt),  
      		INT32_ADDRESS_OFF(patches), gaussPoints, comm_data.nhdgqp, ntriangle, nipp, 
      		comm_data.nearpd, int_data.ws1, ipolator_near);
      	alogger::stopTimer("FMM Initialization");

      	alogger::startTimer("FMM Partitioning");
      	FMM_Partition2(out_size, xb, yb, zb, comm_data.xt, comm_data.yt, comm_data.zt, patches, pointlocs);    
      	alogger::stopTimer("FMM Partitioning");
      	self_metadata meta_data(comm_data.xt, comm_data.yt, comm_data.zt, int_data, vcsio, comm_data.cjvk, ipolymatrix, comm_data.nlqp);
      	alogger::startTimer("Set RHS");
      	setRHS(crhs, out_size, comm_data.xt, comm_data.yt, comm_data.zt, pointlocs);
      	alogger::stopTimer("Set RHS");
      	rj.resize(out_size,0.0);
      	d_complex_t_vec wb(out_size);     
      	for (int i = 0; i < out_size; ++i) {
      		wb[i] = ws[pointlocs[i]] * 0.5 / (4.0 * M_PI);
      	}
      	my_gmres(out_size, ntriangle, nipp, comm_data.nitermax, comm_data.precis, crhs, rj, wb, meta_data, comm_data.fmmVerbose, comm_data.checkFMMDirect, comm_data.gmresRestart, patches, pointlocs, zzsparse, comm_data.mpirank);      
      	if(comm_data.writeTimingOutputs) {
      		printVecToFile(rj, rj.size(), "", create_process_file_name("c_solution_rj",comm_data.mpirank));
      	}
      	return rj;
#else 
      	std::cout << "Dense version is not supported in Partitioned Mode " << std::endl;
      	exit(0); 
#endif
      }
#else
      d_complex_t_vec computeSourceField(int& out_size, int32_vec& patches, int16_vec& pointlocs) {
      	int32_t size = comm_data.init_part_size*nipp;
      	d_complex_t_vec crhs;
      	out_size = size;      
      	if(size > 0){
      		patches.resize(size);
      		pointlocs.resize(size);
      	}
      	for (int i = 0; i < comm_data.init_part_size; ++i) {        
      		for (int j = 0; j < nipp; ++j) {
      			patches[i*nipp + j] = i;
      			pointlocs[i*nipp + j]  = j;
      		}
      	}
#if USE_FMM
      	double eps2 = 0.0;
      	d_vector const& ws = int_data.ws;
      	d_vector_2d vcsio(nipp,d_vector(2,0.0));
      	for (int k = 0; k < nipp; ++k) {
      		vcsio[k][0] = alphao[k] * utils::vos1[0] + betao[k] * utils::vos2[0] + gammao[k] * utils::vos3[0];
      		vcsio[k][1] = alphao[k] * utils::vos1[1] + betao[k] * utils::vos2[1] + gammao[k] * utils::vos3[1];
      	}
      	d_vector_2d gaussPoints;
      	d_vector_2d ipolator_near;
      	double ipolymatrix[nipp*nipp]; 
      	double polymatrix[nipp*nipp];
      	utils::initPolyMatrix(polymatrix, nipp, int_data.alphas, int_data.betas, int_data.gammas, ipolymatrix);      

      	self_metadata meta_data(comm_data.nsupan, comm_data.sunod, int_data, vcsio, comm_data.cjvk, ipolymatrix, comm_data.nlqp);
      	alogger::startTimer("Set Near Quad Points");
      	setNearQuadPoints(gaussPoints, ipolator_near, ipolymatrix, ntriangle, nipp);      
      	alogger::stopTimer("Set Near Quad Points");
      	num_threads(comm_data.fmmAttributes.nthreads);
      	alogger::startTimer("FMM Initialization");
      	FMM_Init(eps2, comm_data.fmmAttributes, size, 
      		DOUBLE_ADDRESS_OFF(xb), DOUBLE_ADDRESS_OFF(yb), DOUBLE_ADDRESS_OFF(zb),  
      		INT32_ADDRESS_OFF(patches), gaussPoints, comm_data.nhdgqp, ntriangle, nipp, 
      		comm_data.nearpd, int_data.ws1, ipolator_near);
      	alogger::stopTimer("FMM Initialization");      
      	alogger::startTimer("FMM Partitioning");
      	FMM_Partition(out_size, xb, yb, zb, patches, pointlocs);    
      	alogger::stopTimer("FMM Partitioning");
      	alogger::startTimer("Set RHS");
      	setRHS(crhs, out_size, patches, pointlocs);
      	alogger::stopTimer("Set RHS");
      	rj.resize(out_size,0.0);				
      	d_complex_t_vec wb(out_size);		 
      	for (int i = 0; i < out_size; ++i) {
      		wb[i] = ws[pointlocs[i]] * 0.5 / (4.0 * M_PI);
      	}
      	my_gmres(out_size, ntriangle, nipp, comm_data.nitermax, comm_data.precis, crhs, rj, wb, meta_data, comm_data.fmmVerbose, comm_data.checkFMMDirect, comm_data.gmresRestart, patches, pointlocs, zzsparse, comm_data.mpirank);      
      	if(comm_data.writeTimingOutputs) {
      		printVecToFile(rj, rj.size(), "", createProcessFile("iterative_solution_rj.dat",comm_data.mpirank));
      	}
#else 			
      	d_complex_t_vec ipiv(size);
      	setRHS(crhs, out_size, patches, pointlocs);
      	d_complex_t_vec zzsparse_temp(zzsparse.size());
				for (int i = 0; i < size; ++i) // transpose
      	  for (int j = 0; j < size; ++j)
      	 	  zzsparse_temp[j * size + i] = zzsparse[i * size + j];     
      	int32_t info;
      	int rhs  = 1;
      	lapack_access::zgesv_(&size,&rhs,(d_complex_t_ptr)&zzsparse_temp[0],&size,
      		(d_complex_t_ptr)&ipiv[0],   (d_complex_t_ptr)&crhs[0],&size,&info);
      	if(info>0) std::cout<<"Error in LU inversion: " << info << std::endl;
      	if(info<0) std::cout<<"Error in LU inputs: "    <<-info << std::endl;
      	printVecToFile(crhs, crhs.size(), "", createProcessFile("lapack_solution_rj.dat",comm_data.mpirank));
      	rj = crhs;
#endif
      	return rj;

      }
#endif
    };
} // end namespace bemfmm
#endif
