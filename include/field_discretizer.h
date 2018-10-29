#ifndef DISCRETIZER
#define DISCRETIZER
#include "utils.h"


namespace bemfmm {
  class field_discretizer {

public:
    static void discretize(integral_data const& int_data, common_data const& comm_data,
                           d_complex_t_vec& zzsparse, d_complex_t_vec& self_contributions) {
      int const& nipp      = comm_data.nipp; 
      int const& ntriangle = comm_data.ntriangle;
      int const& nnpt      = comm_data.nnpt;
      double nearpd        = comm_data.nearpd;
      int** const& nsupan  = comm_data.nsupan;
      double** const& sunod = comm_data.sunod;
      zzsparse.resize((nipp*ntriangle) * (nipp*ntriangle));          
      int_data_type const& alphas = int_data.alphas;
      int_data_type const& betas  = int_data.betas;
      int_data_type const& gammas = int_data.gammas;
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
      d_complex_t    const& cjvk   = comm_data.cjvk;
      int_data_type const& alphao = int_data.alphao;
      int_data_type const& betao  = int_data.betao;
      int_data_type const& gammao = int_data.gammao;
      int_data_type const& alphas1 = int_data.alphas1;
      int_data_type const& betas1  = int_data.betas1;
      int_data_type const& gammas1 = int_data.gammas1;
      int_data_type const& ws1     = int_data.ws1;
      int_data_type const& ws      = int_data.ws;
      int const& nhdgqp = comm_data.nhdgqp;
      int const& nlqp   = comm_data.nlqp;
      d_vector const& wwo = int_data.wwo;
      d_vector const& xxo = int_data.xxo;
      d_vector const& wws = int_data.wws;
      d_vector const& xxs = int_data.xxs;
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
#pragma omp parallel for default(none) shared(zzsparse,self_contributions, xxo, vus, ipolymatrix, ntriangle, i, nipp, nnpt, sunod, nsupan, nearpd, cjvk, I, ws, nhdgqp, ipolator_near, ws1, utils::vos1, utils::vos2, utils::vos3, nlqp, vmmm, xxs, wws, wwo) firstprivate(nodes_sn, vippo, nodeo_sn, vcsio, vipps, vcsis, vvs, vrrr, rrr, ielement00, z00, vcsik_, vgaussk, aa1, aa2, bb1, vvo, vuo, ruv1, ruv2, veuv1, veuv2, jacob_star, ppp, qqq, vcsik, vrr0, rr0, ipolator, polyvector, valf)
        for (int j = 0; j < ntriangle; ++j) {
          nodes_sn = nsupan[j];
          if(i!=j) {        
            for (int k = 0; k < nipp; ++k) {
#if USE_MMAP||USE_PART
#else            
              utils::getR(vcsio[k], nodeo_sn, sunod,  nnpt, vippo);
#endif

              for (int l = 0; l < nipp; ++l) {
#if USE_MMAP||USE_PART
#else
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
#if USE_MMAP||USE_PART
#else
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
#if USE_MMAP||USE_PART
#else
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
  };
} // end namespace bemfmm
#endif
