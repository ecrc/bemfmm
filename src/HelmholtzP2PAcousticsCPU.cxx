//#include "fmm/kernel.h"
#include "fmm/simdvec.h"

namespace exafmm {
  namespace kernel {
//     real_t eps2;
//     complex_t wavek;
//     vec3 Xperiodic;
//     int nhdgqp;                                          
//     int nipp;                                           
// #if EXAFMM_NEARF_TREATMENT 
//     std::vector<std::vector<double> > ipolator_near;     
//     std::vector<real_t> ws;       
//     real_t nearpd;
// #endif    
    void P2P(C_iter Ci, C_iter Cj, bool mutual) {
      real_t wave_r = std::real(wavek);
      real_t wave_i = std::imag(wavek);
      B_iter Bi = Ci->BODY;
      B_iter Bj = Cj->BODY;
      int ni = Ci->NBODY;
      int nj = Cj->NBODY;
      int i = 0;            
#if EXAFMM_USE_SIMD                     
      assert(false); 
      const int isize = ceil((double)ni/BITMAP_SIZE);
      uint16_t locs_bitmap[nj*isize];
      uint16_t ilocs[NSIMD];  
      memset(ilocs, 0, NSIMD * sizeof(uint16_t));
      memset(locs_bitmap, 0, nj * isize * sizeof(uint16_t));
      simdvec wave_rvec = wave_r;
      simdvec wave_ivec = wave_i;
      simdvec zero = 0.0;      
      simdvec one = 1.0;
      for ( ; i<=ni-NSIMD; i+=NSIMD) {
        simdvec tnear = nearpd;
        ksimdvec pot_r = zero;
        ksimdvec pot_i = zero;
        ksimdvec ax_r = zero;
        ksimdvec ax_i = zero;
        ksimdvec ay_r = zero;
        ksimdvec ay_i = zero;
        ksimdvec az_r = zero;
        ksimdvec az_i = zero;

        simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(Bi,i);
        simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(Bi,i);
        simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(Bi,i);
        simdvec mi_r = SIMD<simdvec,4,NSIMD>::setBody(Bi,i);
        simdvec mi_i = SIMD<simdvec,5,NSIMD>::setBody(Bi,i);
        simdvec patchi = SIMD<simdvec,6,NSIMD>::setBody(Bi,i);
        simdvec idx = SIMD<simdvec,0,NSIMD>::setIndex(i+1);

        simdvec dx = Xperiodic[0];
        xi -= dx;
        simdvec dy = Xperiodic[1];
        yi -= dy;
        simdvec dz = Xperiodic[2];
        zi -= dz;

        for (int j=0; j<nj; j++) {
          dx = Bj[j].X[0];
          dx -= xi;
          dy = Bj[j].X[1];
          dy -= yi;
          dz = Bj[j].X[2];
          dz -= zi;
          simdvec patchj = Bj[j].PATCH;
          simdvec R2 = eps2;
          R2 += dx * dx;
          simdvec mj_r = std::real(Bj[j].SRC * Bj[j].QWEIGHT);
          R2 += dy * dy;
          simdvec mj_i = std::imag(Bj[j].SRC * Bj[j].QWEIGHT);
          R2 += dz * dz;
          simdvec invR = rsqrt(R2);
          simdvec R = one / invR;          
          invR &= R2 > zero;
          R &= R2 > zero;
          dx = idx;
          idx &= R < tnear;                            
          idx &= patchj != patchi; 
          for (int k=0; k<NSIMD; k++) { 
            if(idx[k] > 0) {
              locs_bitmap[j*isize+(i+k)/BITMAP_SIZE] ^= 0x1 << (i+k)%BITMAP_SIZE;              
            }
          }
          idx = dx;
          invR &= R > tnear;                          
          R    &= R > tnear;          
          R &= patchj != patchi;     
          invR &= patchj != patchi;
          simdvec tmp = mi_r * mj_r - mi_i * mj_i;
          mj_i = mi_r * mj_i + mi_i * mj_r;
          mj_r = tmp;
          tmp = invR / exp(wave_ivec * R);
          simdvec coef_r = cos(wave_rvec * R) * tmp;
          simdvec coef_i = sin(wave_rvec * R) * tmp;
          tmp = mj_r * coef_r - mj_i * coef_i;
          coef_i = mj_r * coef_i + mj_i * coef_r;
          coef_r = tmp;
          pot_r += coef_r;
          pot_i += coef_i;
        }
        for (int k=0; k<NSIMD; k++) {
          Bi[i+k].TRG[0] += transpose(pot_r, pot_i, k);
        }
      }      
      for (int j = 0; j < nj; ++j) {         
        int simd_count = 0;
        int loc_index = 0;
        real_t mj_r = std::real(Bj[j].SRC);
        real_t mj_i = std::imag(Bj[j].SRC);         
        std::vector<double>& ip_near = ipolator_near[Bj[j].POINT_LOC];
        const int stride = j*isize;
        for (int s = 0; s < isize; ++s) {        
          int element = locs_bitmap[stride+s];
          for (int k = 0; k < BITMAP_SIZE; ++k) {
            if((element & 0x1) == 0x1)
              ilocs[simd_count++] = loc_index; 
            loc_index++;
            element >>= 1;            
            if(simd_count == NSIMD) {
              simdvec xi = SIMD<simdvec,0,NSIMD>::setIndexedBody(Bi,ilocs);
              simdvec yi = SIMD<simdvec,1,NSIMD>::setIndexedBody(Bi,ilocs);
              simdvec zi = SIMD<simdvec,2,NSIMD>::setIndexedBody(Bi,ilocs);
              simdvec mi_r = SIMD<simdvec,4,NSIMD>::setIndexedBody(Bi,ilocs);
              simdvec mi_i = SIMD<simdvec,5,NSIMD>::setIndexedBody(Bi,ilocs);
              simdvec patchi = SIMD<simdvec,6,NSIMD>::setIndexedBody(Bi,ilocs);
              simdvec R2 = eps2;
              simdvec dx = Xperiodic[0];
              xi -= dx;
              simdvec dy = Xperiodic[1];
              yi -= dy;
              simdvec dz = Xperiodic[2];
              zi -= dz;     
              simdvec pot_r = 0.0;  
              simdvec pot_i = 0.0;
              for (int ll = 0; ll < nhdgqp; ++ll) {
                dx = Bj[j].GAUSS_NEAR[ll][0];
                dx -= xi;            
                dy = Bj[j].GAUSS_NEAR[ll][1];
                dy -= yi;            
                dz = Bj[j].GAUSS_NEAR[ll][2];
                dz -= zi;   
                R2 += dx * dx;
                R2 += dy * dy;         
                R2 += dz * dz;    
                simdvec mj_near_r = mj_r * ws[ll] * (real_t)ip_near[ll] * NEAR_FIELD_FACTOR;
                simdvec mj_near_i = mj_i * ws[ll] * (real_t)ip_near[ll] * NEAR_FIELD_FACTOR;
                simdvec invR = rsqrt(R2);            
                simdvec R = one / invR;              
                simdvec tmp = mi_r * mj_near_r - mi_i * mj_near_i;
                mj_near_i = mi_r * mj_near_i + mi_i * mj_near_r;
                mj_near_r = tmp;
                tmp = invR / exp(wave_ivec * R);
                simdvec coef_r = cos(wave_rvec * R) * tmp;
                simdvec coef_i = sin(wave_rvec * R) * tmp;            
                tmp = mj_near_r * coef_r - mj_near_i * coef_i;
                coef_i = mj_near_r * coef_i + mj_near_i * coef_r;
                coef_r = tmp;            
                R2 = eps2;
                pot_i += coef_i;
                pot_r += coef_r;
              }
              simd_count = 0;       
              for (int ss = 0; ss < NSIMD; ++ss) {
                Bi[ilocs[ss]].TRG[0] += transpose(pot_r, pot_i, ss);
              }              
            }
          }          
      }
      if(simd_count > 0) {
        for (int s = 0 ; s < simd_count; s++) {
          real_t pot_r = 0.0;
          real_t pot_i = 0.0;
          real_t mi_r = std::real(Bi[ilocs[s]].SRC * Bi[ilocs[s]].QWEIGHT);
          real_t mi_i = std::imag(Bi[ilocs[s]].SRC * Bi[ilocs[s]].QWEIGHT);
          for (int ll = 0; ll < nhdgqp; ++ll) {
            vec3 dX_near = Bi[ilocs[s]].X - Bj[j].GAUSS_NEAR[ll] - Xperiodic;            
            real_t mj_near_r = mj_r * ws[ll] * (real_t)ip_near[ll] * NEAR_FIELD_FACTOR;
            real_t mj_near_i = mj_i * ws[ll] * (real_t)ip_near[ll] * NEAR_FIELD_FACTOR;
            real_t RR = sqrt(norm(dX_near));
            real_t src2_r = mi_r * mj_near_r - mi_i * mj_near_i;
            real_t src2_i = mi_r * mj_near_i + mi_i * mj_near_r;
            real_t expikr_r_= std::exp(wave_i * RR) * RR;     
            real_t expikr_r = std::cos(wave_r * RR) / expikr_r_;
            real_t expikr_i = std::sin(wave_r * RR) / expikr_r_;   
            pot_r += src2_r * expikr_r - src2_i * expikr_i;
            pot_i += src2_r * expikr_i + src2_i * expikr_r;
          }
          Bi[ilocs[s]].TRG[0] += complex_t(pot_r, pot_i);
        }
      }
    }
#endif
#if EXAFMM_USE_AUTOVEC
#pragma simd
#endif
      for ( ; i<ni; i++) {
        real_t pot_r = 0.0;
        real_t pot_i = 0.0;
        real_t mi_r = std::real(Bi[i].SRC * Bi[i].QWEIGHT);
        real_t mi_i = std::imag(Bi[i].SRC * Bi[i].QWEIGHT);
        for (int j=0; j<nj; j++) { 
#if EXAFMM_USE_AUTOVEC          
          real_t dX[3];
          dX[0] = xi[0] - Bj[j].X[0];
          dX[1] = xi[1] - Bj[j].X[1];
          dX[2] = xi[2] - Bj[j].X[2];
          real_t R2 = dX[0]*dX[0] + dX[1]*dX[1] + dX[2]*dX[2] + eps2;
#else          
          vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;          
          real_t R2 = norm(dX) + eps2;   
#endif
          if(Bi[i].PATCH != Bj[j].PATCH) {
            if (R2 != 0) {
              real_t R = sqrt(R2);
#if EXAFMM_NEARF_TREATMENT
              if(R <= nearpd) {
                real_t mj_r = std::real(Bj[j].SRC);
                real_t mj_i = std::imag(Bj[j].SRC);
                std::vector<double>& ip_near = ipolator_near[Bj[j].POINT_LOC];                
                for (int ll = 0; ll < nhdgqp; ++ll) {
#if EXAFMM_USE_AUTOVEC                            
                  real_t dX_near[3];
                  dX_near[0] = xi[0] - Bj[j].GAUSS_NEAR[ll][0];                  
                  dX_near[1] = xi[1] - Bj[j].GAUSS_NEAR[ll][1];
                  dX_near[2] = xi[2] - Bj[j].GAUSS_NEAR[ll][2];
                  real_t RR = dX_near[0]*dX_near[0] + dX_near[1]*dX_near[1] + dX_near[2]*dX_near[2] + eps2;
                  RR = sqrt(RR);
#else
                  vec3 dX_near = Bi[i].X - Bj[j].GAUSS_NEAR[ll] - Xperiodic;         
                  real_t RR = sqrt(norm(dX_near));   
#endif              
                  real_t mj_near_r = mj_r * ws[ll] * (real_t)ip_near[ll] * NEAR_FIELD_FACTOR;
                  real_t mj_near_i = mj_i * ws[ll] * (real_t)ip_near[ll] * NEAR_FIELD_FACTOR;                                    
                  real_t src2_r = mi_r * mj_near_r - mi_i * mj_near_i;
                  real_t src2_i = mi_r * mj_near_i + mi_i * mj_near_r;
                  real_t expikr_r_= std::exp(wave_i * RR) * RR;     
                  real_t expikr_r = std::cos(wave_r * RR) / expikr_r_;
                  real_t expikr_i = std::sin(wave_r * RR) / expikr_r_;   
                  pot_r += src2_r * expikr_r - src2_i * expikr_i;
                  pot_i += src2_r * expikr_i + src2_i * expikr_r;              
                 }              
              } else 
#endif
               {
                real_t mj_r = std::real(Bj[j].SRC * Bj[j].QWEIGHT);
                real_t mj_i = std::imag(Bj[j].SRC * Bj[j].QWEIGHT);
                real_t src2_r = mi_r * mj_r - mi_i * mj_i;
                real_t src2_i = mi_r * mj_i + mi_i * mj_r;      
                real_t expikr_r_= std::exp(wave_i * R) * R;     
                real_t expikr_r = std::cos(wave_r * R) / expikr_r_;
                real_t expikr_i = std::sin(wave_r * R) / expikr_r_;   
                pot_r += src2_r * expikr_r - src2_i * expikr_i;
                pot_i += src2_r * expikr_i + src2_i * expikr_r;
              }
            }
          }
        }
        Bi[i].TRG[0] += complex_t(pot_r, pot_i);
      }
    }

    void P2P(C_iter C) {
      assert(false);
    }
  }
}
