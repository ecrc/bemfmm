#ifndef kernel_h
#define kernel_h
#include <cmath>
#include "types.h"

namespace exafmm {
  namespace kernel {
    const complex_t I(0.,1.);
    real_t eps2;                                                //!< Epslion squared
    complex_t wavek;                                            //!< Helmholtz wave number
    vec3 Xperiodic;                                             //!< Periodic coordinate offset    
    int nhdgqp;                                                 //!< Number of high degree gauss quadrature points 
    int nipp;                                                   //!< Number of integration points    
#if EXAFMM_NEARF_TREATMENT    
    std::vector<std::vector<double> > ipolator_near;            //!< Basis vector for near interactions
    std::vector<real_t> ws;                                     //!< Gauss quadrature integral points
    real_t nearpd;                                              //!< Minimum near patch distance
#endif    
    void setup();                                               //!< Setup phase for kernels
    void P2P(C_iter Ci, C_iter Cj, bool mutual);                //!< P2P kernel between cells Ci and Cj
    void P2P(C_iter C);                                         //!< P2P kernel for cell C
    void P2M(C_iter C);                                         //!< P2M kernel for cell C
    void M2M(C_iter Ci, C_iter C0);                             //!< M2M kernel for one parent cell Ci
    void M2L(C_iter Ci, C_iter Cj, bool mutual);                //!< M2L kernel between cells Ci and Cj
    void L2L(C_iter Ci, C_iter C0);                             //!< L2L kernel for one child cell Ci
    void L2P(C_iter Ci);                                        //!< L2P kernel for cell Ci
  }
}

#include "../src/HelmholtzP2PAcousticsCPU.cxx"
#include "../src/HelmholtzSphericalCPU.cxx"

#endif
