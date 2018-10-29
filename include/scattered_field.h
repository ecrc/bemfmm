#ifndef FAR_NEAR_FIELD
#define FAR_NEAR_FIELD
#include "utils.h"
#include "global_data.h"
namespace bemfmm {    
  enum field_mode
  {
    FAR = 0,
    NEAR_NEAR_SCHEME = 1,
    NEAR_FAR_SCHEME = 2
  };

  class scattered_field {
  private:        
    integral_data const& int_data;
    common_data const& comm_data;            
    d_complex_t_vec const& rj;
    int** nsupan;
    int nipp; 
    int const nnpt;
    int const ntriangle;
    int_data_type const alphas;
    int_data_type const betas;
    int_data_type const gammas;

  public:
    scattered_field(integral_data const& _int_data, common_data const& _comm_data, d_complex_t_vec const& _rj, int** _nsupan):
                  int_data(_int_data),comm_data(_comm_data), rj(_rj), 
                  nsupan(_nsupan), 
                  nipp(comm_data.nipp), nnpt(comm_data.nnpt), ntriangle(comm_data.ntriangle),
                  alphas(int_data.alphas), betas(int_data.betas),
                  gammas(int_data.gammas) { }


    
    void calculateScatteredField(int mode, int mpirank, int mpisize, int partition_size, bool write_output, int32_vec const& patches, int16_vec const& pointlocs, std::string file_name) {     
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
  
  };
}

#endif
