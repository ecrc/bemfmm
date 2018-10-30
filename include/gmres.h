#ifndef GMRES
#define GMRES
#include "utils.h"
#include "fmm/acoustics_wrapper.h"
#include "self_metadata.h"

using namespace bemfmm;
using namespace bemfmm::utils;
void givens_loc(d_complex_t const& z1, d_complex_t const& z2, double& c, d_complex_t& s) {
  double abs_z1 = std::abs(z1);
  double abs_z2 = std::abs(z2);
  double vnormz = sqrt(abs_z1*abs_z1 + abs_z2*abs_z2); 
  if(abs_z1!= 0) {
    c = abs_z1/vnormz;
    s = z2 / z1 * c;
  } else if(abs_z2 != 0) {
    c = 0.0;
    s = z2 / abs_z2;
  } else {
    c = 1.0;
    s = 0.0;
  }
}

double scnrm22(int nd, d_complex_t_vec const& cy9) {
  d_complex_t result = 0.0;
  for (int i = 0; i < nd; ++i) {
    result += std::conj(cy9[i]) * cy9[i];    
  }
  d_complex_t sum; 
  MPI_Allreduce(&result, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  result = sqrt(sum);
  return std::real(result);
}

void cutrisub(int n, d_complex_t_vec_2d const& cma, d_complex_t_vec cvx, d_complex_t_vec& cvy) {
  for (int i = n-1; i >= 0; --i) {
    cvy[i] = cvx[i];
    for (int j = i+1; j < n; ++j) {
      cvy[i]-= cma[j][i] * cvy[j];      
    }    
    cvy[i] /= cma[i][i]; 
  }
}

template <typename VecArr>
void addSelfCorrections(VecArr& target, VecArr const& source, int partition_size, int32_vec const& patches, d_complex_t_vec const& self, int ntriangle, int nipp) { 
  int selfIndex = 0;
  for (int i = 0; i < partition_size; ++i) {    
    for (int j = 0; j < nipp; ++j) {
      target[i] += source[(i/nipp)*nipp + j]*self[selfIndex++];
    }    
  }
}
#if USE_PART
void buildSingularityCorrections(self_metadata const& metadata, int partition_size, int32_vec const& patches, int16_vec const& pointlocs, int ntriangle, int nipp, d_complex_t_vec & self) {
  self.resize(partition_size*nipp, 0);
  double ppp,qqq,aa1[3]={0.0, 0.0, 0.0},aa2[3]={0.0, 0.0, 0.0},bb1[3]={0.0, 0.0, 0.0},vmmm[3]={0.0, 0.0, 0.0},valf[3]={0.0, 0.0, 0.0};
  double veuv1[2]={0.0, 0.0},veuv2[2]={0.0, 0.0};
  double vcsik[2]={0.0, 0.0};
  double vuo[3]={0.0, 0.0, 0.0},vvo[3]={0.0, 0.0, 0.0};  
  d_complex_t I(0.,1.);
  d_vector const& wwo = metadata.wwo;
  d_vector const& xxo = metadata.xxo;
  d_vector const& wws = metadata.wws;
  d_vector const& xxs = metadata.xxs;
  d_vector const& x   = metadata.x;
  d_vector const& y   = metadata.y;
  d_vector const& z   = metadata.z;
  int const& nlqp     = metadata.nlqp;
  double const* ipolymatrix = metadata.ipolymatrix;
  d_vector_2d const& vcsio = metadata.vcsio;
  d_complex_t const& cjvk   = metadata.cjvk;
  double vrrr[3]={0.0, 0.0, 0.0};
  double vrr0[3]={0.0, 0.0, 0.0};  
  d_ptr ruv1, ruv2;
#pragma omp parallel for private(vrrr, vrr0, ruv1, ruv2, ppp, qqq, aa1, aa2, bb1, vmmm, valf, veuv1, veuv2, vcsik, vuo, vvo)
  for (int i = 0; i < partition_size; ++i) {
    const int patch = i/nipp;
    const int index = patch*nipp;
  //  for(int nn = 0; nn < 6; nn++) {
  //    std::cout << x[index+nn] << " , " << y[index+nn] << " , " << z[index+nn] << std::endl ;
  //  }

    aa1[0]=(x[index+0]+x[index+2]-2.0*x[index+5]);
    aa2[0]=(x[index+2]+x[index+3]-x[index+4]-x[index+5]);
    bb1[0]=(x[index+1]+x[index+2]-2.0*x[index+4]);
    
    aa1[1]=(y[index+0]+y[index+2]-2.0*y[index+5]);
    aa2[1]=(y[index+2]+y[index+3]-y[index+4]-y[index+5]);
    bb1[1]=(y[index+1]+y[index+2]-2.0*y[index+4]);
    
    aa1[2]=(z[index+0]+z[index+2]-2.0*z[index+5]);
    aa2[2]=(z[index+2]+z[index+3]-z[index+4]-z[index+5]);
    bb1[2]=(z[index+1]+z[index+2]-2.0*z[index+4]);   
    double polyvector[nipp], rrr, rr0;
    int k = pointlocs[i];
    //for (int k = 0; k < nipp; ++k) 
    {      
      utils::getUV_vector(vcsio[k], x, y, z, patch, nipp, vuo, vvo); 
      for (int l = 0; l < nipp; ++l) {
        d_complex_t z00 = 0.0;
        for (int d = 0; d < 3; ++d) {
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
          double jacob_star = std::abs(veuv1[0]*veuv2[1]-veuv1[1]*veuv2[0]);
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
              double ipolator = 0.0;
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
        self[(i*nipp) + l] += z00;               
      }      
    }
  }
}
#else
void buildSingularityCorrections(self_metadata const& metadata, int partition_size, int32_vec const& patches, int16_vec const& pointlocs, int ntriangle, int nipp, d_complex_t_vec & self) {
  self.resize(partition_size*nipp, 0);
  double ppp,qqq,aa1[3]={0.0, 0.0, 0.0},aa2[3]={0.0, 0.0, 0.0},bb1[3]={0.0, 0.0, 0.0},vmmm[3]={0.0, 0.0, 0.0},valf[3]={0.0, 0.0, 0.0};
  double veuv1[2]={0.0, 0.0},veuv2[2]={0.0, 0.0};
  double vcsik[2]={0.0, 0.0};
  double vuo[3]={0.0, 0.0, 0.0},vvo[3]={0.0, 0.0, 0.0};  
  d_complex_t I(0.,1.);
  d_vector const& wwo = metadata.wwo;
  d_vector const& xxo = metadata.xxo;
  d_vector const& wws = metadata.wws;
  d_vector const& xxs = metadata.xxs;
  int const& nlqp     = metadata.nlqp;
  double const* ipolymatrix = metadata.ipolymatrix;
  int** const& nsupan = metadata.nsupan;
  double** const& sunod   = metadata.sunod;
  d_vector_2d const& vcsio = metadata.vcsio;
  d_complex_t const& cjvk   = metadata.cjvk;
  double vrrr[3]={0.0, 0.0, 0.0};
  double vrr0[3]={0.0, 0.0, 0.0};  
  d_ptr ruv1, ruv2;
#pragma omp parallel for private(vrrr, vrr0, ruv1, ruv2, ppp, qqq, aa1, aa2, bb1, vmmm, valf, veuv1, veuv2, vcsik, vuo, vvo)
  for (int i = 0; i < partition_size; ++i) {
    int* nodeo_sn = nsupan[patches[i]];
    double polyvector[nipp], rrr, rr0;
 
    for (int d = 0; d < 3; ++d) {
      aa1[d]=(sunod[d][nodeo_sn[0]]+sunod[d][nodeo_sn[2]]-2.0*sunod[d][nodeo_sn[5]]);
      aa2[d]=(sunod[d][nodeo_sn[2]]+sunod[d][nodeo_sn[3]]-sunod[d][nodeo_sn[4]]-sunod[d][nodeo_sn[5]]);
      bb1[d]=(sunod[d][nodeo_sn[1]]+sunod[d][nodeo_sn[2]]-2.0*sunod[d][nodeo_sn[4]]);
    }
    int k = pointlocs[i];
    //for (int k = 0; k < nipp; ++k) 
    {      
      utils::getUV_vector(vcsio[k], nodeo_sn, sunod, nipp, vuo, vvo); 
      for (int l = 0; l < nipp; ++l) {
        d_complex_t z00 = 0.0;
        for (int d = 0; d < 3; ++d) {
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
          double jacob_star = std::abs(veuv1[0]*veuv2[1]-veuv1[1]*veuv2[0]);
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
              double ipolator = 0.0;
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
        self[(i*nipp) + l] += z00;               
      }      
    }
  }
}
#endif

int my_gmres(int partition_size, int ntriangle, int nipp, int nitermax, double precis, d_complex_t_vec& crhs, d_complex_t_vec& rj, 
 d_complex_t_vec& weights, self_metadata const& metadata, int fmm_verbose, bool verify, int gmresRestart, int32_vec const& patches, int16_vec const& pointlocs,
 d_complex_t_vec const& zzsparse,int mpirank) {

const int m = gmresRestart;
const int nd = partition_size;
//const int sample_direct = 300;
const int sample_direct = 500;
int itertotal, iterout, n, k;
d_complex_t_vec& cx = rj;
d_complex_t_vec cb(m+1,0.0);
d_complex_t_vec chr(m+1,0.0);  
d_complex_t_vec self;
d_complex_t_vec_2d cy(m+1, d_complex_t_vec(nd,0.0));
d_complex_t_vec_2d ch(m,   d_complex_t_vec(m+1,0.0));
d_complex_t_vec_2d cy_t(nd,d_complex_t_vec(m,0.0));
d_complex_t_vec test(nd ,0);  
std::vector<int> sample_addresses(sample_direct ,0);  
buildSingularityCorrections(metadata, nd, patches, pointlocs, ntriangle, nipp, self);
bool init = true;
d_complex_t_vec cs(m,0.0);
d_vector rc(m,0.0);
d_complex_t ctemp = 0.0;
double ay0, be, bea;
alogger::startTimer("GMRES Time");
ay0   = scnrm22(nd, crhs);
itertotal = 0;
for (iterout = 0; iterout < nitermax; ++iterout) {
  itertotal++; 
  if(!init) {
    alogger::startTimer("FMM Time");
    exafmm::FMM_B2B(ZDOUBLE_ADDRESS_OFF(cy[m]), ZDOUBLE_ADDRESS_OFF(cx), ZDOUBLE_ADDRESS_OFF(weights), fmm_verbose);      
    addSelfCorrections(cy[m], cx, nd, patches, self, ntriangle, nipp);
    alogger::stopResetTimer("FMM Time");        
  }
  else init = false;  
  for (int i = 0; i < nd; ++i) cy[m][i]  = crhs[i] - cy[m][i];
    cb[0] = scnrm22(nd,cy[m]);
  be = std::real(cb[0]/ay0);
  if(alogger::verbose) {
		std::cout << std::setw(alogger::stringLength) << std::fixed << std::left
    << "Iter Outer" << " : " << iterout << std::endl;			
		std::cout << std::setw(alogger::stringLength) << std::fixed << std::left
    << "Iter Total" << " : " << itertotal << std::endl;
		std::cout << std::setw(alogger::stringLength) << std::fixed << std::left
    << "GMRES Residual Norm" << " : " << std::setprecision(9) << std::scientific<< be <<std::endl;      
		std::cout << "==========================================================" << std::fixed << std::left <<std::endl;		 		
  }
  if(be < precis || itertotal > nitermax) break;
  for (int i = 0; i < nd; ++i) cy[0][i] = cy[m][i]/cb[0];
    for (n = 0; n < m; ++n) {
      alogger::startTimer("Iter Time");
      itertotal++;       
      alogger::startTimer("FMM Time");
      exafmm::FMM_B2B(ZDOUBLE_ADDRESS_OFF(cy[n+1]), ZDOUBLE_ADDRESS_OFF(cy[n]), ZDOUBLE_ADDRESS_OFF(weights), fmm_verbose);                               
      if(verify) {
				alogger::stopTimer("Solving AX=B",0);
				exafmm::DirectSample(sample_direct, ZDOUBLE_ADDRESS_OFF(test), sample_addresses);												 
				alogger::logNumericalError(test, cy[n+1], sample_addresses, sample_direct ,"FMM vs. Direct");
				alogger::startTimer("Solving AX=B");
			}
      addSelfCorrections(cy[n+1], cy[n], partition_size, patches, self, ntriangle, nipp);        
      alogger::stopResetTimer("FMM Time");        
      for (k = 0; k <= n; ++k) {
        dot_product(cy[k],cy[n+1],ch[n][k],nd);                             
      }
      MPI_Allreduce(ZDOUBLE_ADDRESS_OFF(ch[n]), ZDOUBLE_ADDRESS_OFF(chr), n+1,
        MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
      for (k = 0; k <= n; ++k) { 
        ch[n][k] = chr[k];
        for (int i = 0; i < nd; ++i) cy[n+1][i]-=cy[k][i]*ch[n][k];
      }  

    ch[n][n+1] = scnrm22(nd,cy[n+1]);
    if(n < m) for (int i = 0; i < nd; ++i) cy[n+1][i]/=ch[n][n+1];
    if(n != 0) {        
      for (int k = 0; k <= n-1; ++k) {          
        ctemp    = rc[k]*ch[n][k+1] - cs[k]*ch[n][k];
        ch[n][k] = rc[k]*ch[n][k]   + std::conj(cs[k]) * ch[n][k+1];
        ch[n][k+1]   = ctemp;
      }
    }      
    givens_loc(ch[n][n],ch[n][n+1],rc[n],cs[n]);
    cb[n+1] =  -cs[n]*cb[n];
    cb[n]   =   rc[n]*cb[n];  
    ch[n][n]   = rc[n]*ch[n][n]+ std::conj(cs[n])*ch[n][n+1];
    ch[n][n+1] = 0;
    bea = std::abs(cb[n+1])/ay0;
    if(alogger::verbose) {
		  std::cout << std::setw(alogger::stringLength) << std::fixed << std::left
      << "Iteration" << " : " << itertotal << std::endl;
			std::cout << std::setw(alogger::stringLength) << std::fixed << std::left
      << "GMRES Residual Norm" << " : " << std::setprecision(9) << std::scientific<< bea <<std::endl;      
			std::cout << "==========================================================" << std::fixed << std::left <<std::endl;     
		}    
    if(n == m-1 || bea < precis) {          
      cutrisub(n+1, ch, cb, cb);     
      transpose(cy,m, nd, cy_t);     
      matVecMul(cy_t,cb,nd,n+1,cy[m]);       
      for (int i = 0; i < nd; ++i) cx[i] = cx[i] + cy[m][i];  
        break;     
    }
    alogger::stopResetTimer("Iter Time");
  }
alogger::stopTimer("GMRES Time");
}
return 0;
}


#endif
