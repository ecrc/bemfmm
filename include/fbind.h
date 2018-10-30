/**
 *
 * @file fbind.h
 *
 * @copyright 2018 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
   * @file fbind.h   
   *
   *  BEMFMM is a software package provided by King Abdullah University of Science and Technology (KAUST)
   *
   * @version 1.0.0
   * @author Mustafa Abduljabbar & Mohammed Al Farhan
   *
 **/

#ifndef FBIND_H
#define FBIND_H
#include <complex>
namespace bemfmm {
    namespace lapack_access {
#ifdef __cplusplus
  extern "C" {    
    void zgesv_(int *lengthA,    int *widthF, std::complex<double> *A,      int *leadingDemA, 
                std::complex<double> *permMat, std::complex<double> *B, int *leadingDemB, int *errorCheck);
    void dgesv_(int *lengthA,    int *widthF, double *A,      int *leadingDemA, 
                double *permMat, double *B, int *leadingDemB, int *errorCheck); 
  }
#endif
#ifdef __cplusplus
#endif
    }
}
#endif
