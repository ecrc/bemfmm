/**
 *
 * @file utils.h
 *
 * @copyright 2018 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 *
 * @author Mustafa Abduljabbar [mustafa.abduljabbar@kaust.edu.sa] and Mohammed Al Farhan [mohammed.farhan@kaust.edu.sa].
 *
 **/

#ifndef UTILS
#define UTILS
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <limits>
#include <math.h>
#include <vector>
#include <fstream>
#include "fbind.h"
#include "logger.h"
#include "args.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "thread.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <assert.h>
#include "fmm/logger.h"
using namespace exafmm;
#define ROUND_OFF(X,DIG) round(X*DIG)/DIG

namespace bemfmm {
  const int dim = 3;
  const int FILE_NAME_SIZE = 40;
  const int MASTER_NODE_ID = 0;
  const size_t INT_BYTE_LIMIT = (((size_t) 2) << ((((sizeof(int) * 8)))-5));
  typedef std::complex<double> d_complex_t;
  typedef d_complex_t* d_complex_t_ptr;
  typedef std::vector<d_complex_t> d_complex_t_vec;
  typedef std::vector<d_complex_t_vec> d_complex_t_vec_2d;
  typedef std::vector<double> d_vector;     
  typedef std::vector<d_vector> d_vector_2d;    
  typedef d_vector int_data_type; 
  typedef std::vector<int32_t> int32_vec;
  typedef std::vector<int16_t> int16_vec;
  typedef std::vector<int32_vec> int32_vec_2d;
  typedef double* d_ptr;
  typedef uint32_t* uint_ptr;
  typedef int32_t* int_ptr;  
  int node_id;  
  struct fmm_context {
    d_complex_t_vec const& weights;
    d_complex_t_vec const& self_corrections;
    int32_vec const& patches;
    int partition_size;
    bool fmm_verbose;
    int ntriangle; 
    int nipp;
    bool init;
    public:
      fmm_context(d_complex_t_vec const& w, d_complex_t_vec const& self, int32_vec const& pats, int partition, bool verbose, int ntri, int pp): 
                weights(w), self_corrections(self), patches(pats), partition_size(partition), fmm_verbose(verbose), 
                ntriangle(ntri), nipp(pp), init(false) {} 
  };

  namespace utils {
    double vos1[2] = {1.0, 0.0};
    double vos2[2] = {0.0, 1.0};
    double vos3[2] = {0.0, 0.0};
    
#if USE_PART
    template<typename ArrT1, typename ArrT2, typename OutArr>
    inline void getR(ArrT1 const& vkesi, ArrT2 const& sunod, int const& nnpt, OutArr& rr, int dim = 3) {
      double vkesi3 = 1.0 - vkesi[0] - vkesi[1];
      if(nnpt == 6) {
        for (int i = 0; i < dim; ++i) {
          rr[i] = vkesi[0]*(2.0*vkesi[0]-1.0)*sunod[0*3+i]
             +vkesi[1]*(2.0*vkesi[1]-1.0)*sunod[1*3+i]
             +vkesi3*(2.0*vkesi3-1.0)*sunod[2*3+i]
             +4.0*vkesi[0]*vkesi[1]*sunod[3*3+i]
             +4.0*vkesi[1]*vkesi3*sunod[4*3+i]
             +4.0*vkesi3*vkesi[0]*sunod[5*3+i];
        }
      } else {
        std::cout<<"can not get r for other order mesh now"<<std::endl;
      }
    }

    template<typename ArrT1, typename ArrT2, typename OutArr>
    inline void getR(ArrT1 const& vkesi, ArrT2 const& x, ArrT2 const& y, ArrT2 const& z, int const& patch, int const& nnpt, OutArr& rr, int dim = 3) {
      double vkesi3 = 1.0 - vkesi[0] - vkesi[1];
      if(nnpt == 6) {        
        const int index = patch*nnpt;
        if((index + 5) > x.size()) {
          std::cout << index + 5 << "  " << x.size() << "  " << std::endl;
          assert(false);
        }
        rr[0] = vkesi[0]*(2.0*vkesi[0]-1.0)*x[index+0]
           +vkesi[1]*(2.0*vkesi[1]-1.0)*x[index+1]
           +vkesi3*(2.0*vkesi3-1.0)*x[index+2]
           +4.0*vkesi[0]*vkesi[1]*x[index+3]
           +4.0*vkesi[1]*vkesi3*x[index+4]
           +4.0*vkesi3*vkesi[0]*x[index+5];

        rr[1] = vkesi[0]*(2.0*vkesi[0]-1.0)*y[index+0]
           +vkesi[1]*(2.0*vkesi[1]-1.0)*y[index+1]
           +vkesi3*(2.0*vkesi3-1.0)*y[index+2]
           +4.0*vkesi[0]*vkesi[1]*y[index+3]
           +4.0*vkesi[1]*vkesi3*y[index+4]
           +4.0*vkesi3*vkesi[0]*y[index+5];
           
        rr[2] = vkesi[0]*(2.0*vkesi[0]-1.0)*z[index+0]
           +vkesi[1]*(2.0*vkesi[1]-1.0)*z[index+1]
           +vkesi3*(2.0*vkesi3-1.0)*z[index+2]
           +4.0*vkesi[0]*vkesi[1]*z[index+3]
           +4.0*vkesi[1]*vkesi3*z[index+4]
           +4.0*vkesi3*vkesi[0]*z[index+5]; 
      } else {
        std::cout<<"can not get r for other order mesh now"<<std::endl;
      }
    }
    
    template<typename ArrT1, typename ArrT2, typename OutArr>
    inline void getUV_vector(ArrT1 const& vkesi, ArrT2 const& x, ArrT2 const& y, ArrT2 const& z, int const& patch, int const& nnpt, OutArr& vi1, OutArr& vi2) {
      if(nnpt == 6) {
        const int index = patch*nnpt;
        vi1[0]=4.0*x[index+5]-3.0*x[index+2]-x[index+0]    //vector u
             +4.0*(x[index+2]+x[index+3]-x[index+4]-x[index+5])*vkesi[1]
             +4.0*(x[index+0]+x[index+2]-2.0*x[index+5])*vkesi[0];
        vi2[0]=4.0*x[index+4]-3.0*x[index+2]-x[index+1]    //vector v
             +4.0*(x[index+2]+x[index+3]-x[index+4]-x[index+5])*vkesi[0]
             +4.0*(x[index+1]+x[index+2]-2.0*x[index+4])*vkesi[1];

        vi1[1]=4.0*y[index+5]-3.0*y[index+2]-y[index+0]    //vector u
             +4.0*(y[index+2]+y[index+3]-y[index+4]-y[index+5])*vkesi[1]
             +4.0*(y[index+0]+y[index+2]-2.0*y[index+5])*vkesi[0];
        vi2[1]=4.0*y[index+4]-3.0*y[index+2]-y[index+1]    //vector v
             +4.0*(y[index+2]+y[index+3]-y[index+4]-y[index+5])*vkesi[0]
             +4.0*(y[index+1]+y[index+2]-2.0*y[index+4])*vkesi[1];

        vi1[2]=4.0*z[index+5]-3.0*z[index+2]-z[index+0]    //vector u
             +4.0*(z[index+2]+z[index+3]-z[index+4]-z[index+5])*vkesi[1]
             +4.0*(z[index+0]+z[index+2]-2.0*z[index+5])*vkesi[0];
        vi2[2]=4.0*z[index+4]-3.0*z[index+2]-z[index+1]    //vector v
             +4.0*(z[index+2]+z[index+3]-z[index+4]-z[index+5])*vkesi[0]
             +4.0*(z[index+1]+z[index+2]-2.0*z[index+4])*vkesi[1];
      } else {
        std::cout<<"can not get r for other order mesh now"<<std::endl;
      }
    }

    template<typename T>
    void fetchCoords(int ntriangle, int nipp, int nnpt, 
                        int_data_type const alphas, int_data_type const betas, 
                        int_data_type const gammas, T & comm_data,
                        d_vector& x, d_vector& y, d_vector& z) {      
      double vos1[2] = {1.0,0.0};
      double vos2[2] = {0.0,1.0};
      double vos3[2] = {0.0,0.0};
      double vcsis[2], vipps[3];
      size_t index = 0;
      for (int i = comm_data.start_triangle; i < comm_data.end_triangle; ++i) {
        double* _sunod = comm_data.triangles + (i*nipp*3);
        for (int j = 0; j < nipp; ++j) {
           vcsis[0] = alphas[j]*vos1[0] + betas[j]*vos2[0] + gammas[j]*vos3[0];
           vcsis[1] = alphas[j]*vos1[1] + betas[j]*vos2[1] + gammas[j]*vos3[1];
           getR(vcsis, _sunod, nnpt, vipps);
           x[index] = vipps[0];
           y[index] = vipps[1];
           z[index++] = vipps[2];           
        }
      }
    }
#else
    template<typename ArrT1, typename ArrT2, typename ArrT3, typename OutArr>
    inline void getR(ArrT1 const& vkesi, ArrT2 const& array, ArrT3 const& sunod, int const& nnpt, OutArr& rr, int dim = 3) {
      double vkesi3 = 1.0 - vkesi[0] - vkesi[1];
      if(nnpt == 6) {
        for (int i = 0; i < dim; ++i) {
          rr[i] = vkesi[0]*(2.0*vkesi[0]-1.0)*sunod[i][array[0]]
             +vkesi[1]*(2.0*vkesi[1]-1.0)*sunod[i][array[1]]
             +vkesi3*(2.0*vkesi3-1.0)*sunod[i][array[2]]
             +4.0*vkesi[0]*vkesi[1]*sunod[i][array[3]]
             +4.0*vkesi[1]*vkesi3*sunod[i][array[4]]
             +4.0*vkesi3*vkesi[0]*sunod[i][array[5]];
        }

      } else {
        std::cout<<"can not get r for other order mesh now"<<std::endl;
      }
    }
    
    template<typename Arr>
    void fetchCoords(int ntriangle, int nipp, int nnpt, int** const& _nsupan, 
                        int_data_type const alphas, int_data_type const betas, 
                        int_data_type const gammas, Arr const& sunod,
                        d_vector& x, d_vector& y, d_vector& z) {      
      double vos1[2] = {1.0,0.0};
      double vos2[2] = {0.0,1.0};
      double vos3[2] = {0.0,0.0};
      int* nodes_sn;
      double vcsis[2], vipps[3];
      size_t index = 0;
      for (int i = 0; i < ntriangle; ++i) {
        nodes_sn = _nsupan[i];
        for (int j = 0; j < nipp; ++j) {
           vcsis[0] = alphas[j]*vos1[0] + betas[j]*vos2[0] + gammas[j]*vos3[0];
           vcsis[1] = alphas[j]*vos1[1] + betas[j]*vos2[1] + gammas[j]*vos3[1];
           getR(vcsis, nodes_sn, sunod, nnpt, vipps);
           x[index] = vipps[0];
           y[index] = vipps[1];
           z[index++] = vipps[2];           
        }
      }
    }

#endif    
    bool endsWith(const std::string& str, const std::string& suffix) {
      return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
    }
 

    size_t getFileSize(const char* filename) {
        struct stat st;
        stat(filename, &st);
        return st.st_size;
    }

		inline void initBasisVector(int ngprcs, int nipp, int_data_type _alphas, 
                                          int_data_type _betas, int_data_type _gammas, 
                                          double* ipolymatrix , d_vector_2d& ipolator,d_vector_2d& vcsik) {      
      d_vector_2d polyvector(ngprcs, d_vector(nipp,0));
      for (int m = 0; m < ngprcs; ++m) {
        vcsik[m][0] = _alphas[m]*utils::vos1[0]+_betas[m]*utils::vos2[0]+_gammas[m]*utils::vos3[0];
        vcsik[m][1] = _alphas[m]*utils::vos1[1]+_betas[m]*utils::vos2[1]+_gammas[m]*utils::vos3[1];
        if(nipp == 6) {
          polyvector[m][0]=1.0;
          polyvector[m][1]=vcsik[m][0];
          polyvector[m][2]=vcsik[m][1];
          polyvector[m][3]=vcsik[m][0]*vcsik[m][1];
          polyvector[m][4]=vcsik[m][0]*vcsik[m][0];
          polyvector[m][5]=vcsik[m][1]*vcsik[m][1];
        }
      }
      for (int i = 0; i < nipp; ++i) {
        for (int m = 0; m < ngprcs; ++m) {
          if(nipp == 6) {
            ipolator[i][m] = 0.0;
            for (int kkk = 0; kkk < nipp; ++kkk) {
              ipolator[i][m] += ipolymatrix[kkk*nipp+i]*polyvector[m][kkk];
            }
          }
        }
      }
    }

    inline void initPolyMatrix(double* const& polymatrix, int const& nipp,
                                       int_data_type const& alphas, int_data_type const& betas, 
                                       int_data_type const& gammas, double* ipolymatrix) {
      double* polymatrix_temp = new double[nipp*nipp];
      double* ipiv = new double[nipp];
      int info;
      d_vector vcsis(2);
      if(nipp == 6) {
        for (int i = 0; i < nipp; ++i) {
          vcsis[0] = alphas[i]*utils::vos1[0] + betas[i]*utils::vos2[0] + gammas[i]*utils::vos3[0];
          vcsis[1] = alphas[i]*utils::vos1[1] + betas[i]*utils::vos2[1] + gammas[i]*utils::vos3[1];
          const uint32_t offset = i*nipp;
          polymatrix_temp[offset+0]=polymatrix[offset+0]=1.0;
          polymatrix_temp[offset+1]=polymatrix[offset+1]=vcsis[0];
          polymatrix_temp[offset+2]=polymatrix[offset+2]=vcsis[1];
          polymatrix_temp[offset+3]=polymatrix[offset+3]=vcsis[0]*vcsis[1];
          polymatrix_temp[offset+4]=polymatrix[offset+4]=vcsis[0]*vcsis[0];
          polymatrix_temp[offset+5]=polymatrix[offset+5]=vcsis[1]*vcsis[1];
        }        
        for (int i = 0; i < nipp; ++i) {
          for (int j = 0; j < nipp; ++j) {
            if(i==j)ipolymatrix[i*nipp+j] = 1.0;
            else ipolymatrix[i*nipp+j] = 0.0;
          }          
        }
        int _nipp = nipp;
        lapack_access::dgesv_(&_nipp, &_nipp,polymatrix_temp,&_nipp,ipiv,ipolymatrix,&_nipp,&info);        
        if(info < 0) {
          std::cerr<< " argument: " << -info << " had illegal value " << std::endl;
        }    
        else if(info > 0){
          std::cerr<< "U(" << info << "," << info<<") is exactly zero"<<std::endl;
        }
      }
      delete[] polymatrix_temp;
      delete[] ipiv;
    }
 
    template<typename ArrT1, typename ArrT2, typename ArrT3, typename OutArr>
    inline void getUV_vector(ArrT1 const& vkesi, ArrT2 const& array, ArrT3 const& sunod, int const& nnpt, OutArr& vi1, OutArr& vi2, int dim = 3) {
      if(nnpt == 6) {
        for (int i = 0; i < dim; ++i) {
          vi1[i]=4.0*sunod[i][array[5]]-3.0*sunod[i][array[2]]-sunod[i][array[0]]    //vector u
            +4.0*(sunod[i][array[2]]+sunod[i][array[3]]-sunod[i][array[4]]-sunod[i][array[5]])*vkesi[1]
            +4.0*(sunod[i][array[0]]+sunod[i][array[2]]-2.0*sunod[i][array[5]])*vkesi[0];   

          vi2[i]=4.0*sunod[i][array[4]]-3.0*sunod[i][array[2]]-sunod[i][array[1]]    //vector v
            +4.0*(sunod[i][array[2]]+sunod[i][array[3]]-sunod[i][array[4]]-sunod[i][array[5]])*vkesi[0]
            +4.0*(sunod[i][array[1]]+sunod[i][array[2]]-2.0*sunod[i][array[4]])*vkesi[1]; 
        }

      } else {
        std::cout<<"can not get r for other order mesh now"<<std::endl;
      }
    }
    template<typename ArrT>
    inline void cross_product(ArrT const& va, ArrT const& vb, ArrT& vc) {
      vc[0] = va[1] * vb[2] - va[2] * vb[1]; 
      vc[1] = va[2] * vb[0] - va[0] * vb[2];
      vc[2] = va[0] * vb[1] - va[1] * vb[0];   
    }

    template<typename Arr1, typename Arr2, typename T>
    inline void dot_product(Arr1 const& va, Arr2 const& vb, T& vc, int dim = 3) {
      vc = va[0] * vb[0];
      for (int i = 1; i < dim; ++i)
        vc += va[i] * vb[i];
    }

    inline void dot_product(d_complex_t_vec const& va, d_complex_t_vec const& vb, d_complex_t& vc, int dim = 3) {
      vc = std::conj(va[0]) * vb[0];
      for (int i = 1; i < dim; ++i)
        vc += std::conj(va[i]) * vb[i];
    }

    template<typename MatT, typename VecT> 
    inline void matVecMul(MatT const& mat, VecT const& vec, int rows, int cols, VecT& res){
#pragma omp parallel for     
      for (int i = 0; i < rows; ++i) {        
        res[i] = mat[i][0]*vec[0];
        for (int j = 1; j < cols; ++j) {
          res[i] += mat[i][j]*vec[j];
        }      
      }
    }

    template<typename MatT> 
    inline void transpose(MatT const& mat, int rows, int cols, MatT& o_mat) {
#pragma omp parallel for      
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
          o_mat[j][i] = mat[i][j];
        }
      }       
    }

    template<typename ArrT, typename T>
    inline void rootDotProduct(ArrT const& va, ArrT const& vb, T& vc, int dim = 3) {
      vc = va[0] * vb[0];
      for (int i = 1; i < dim; ++i)
        vc += va[i] * vb[i];
      vc = sqrt(vc);
    }

    template<typename ArrT, typename T>
    inline void getJacobian(ArrT const& vi1, ArrT const& vi2, T& vii, int dim = 3) { 
      T g11; dot_product(vi1, vi1, g11);
      T g22; dot_product(vi2, vi2, g22);
      T g12; dot_product(vi2, vi1, g12);
      vii = sqrt(g11 * g22 - g12 * g12);
    }
 

    
    template<typename Vec> 
    inline void printVecToFile(Vec const& v, size_t size, const char* headings, std::string file_name) {
      std::fstream fs;
      fs.open(file_name.c_str(), std::fstream::out);
      //fs<<headings<<std::endl;
      for (size_t i = 0; i < size; ++i) {
        fs << std::setprecision(alogger::precis) << std::real(v[i]) << " " 
                      << std::imag(v[i]) << std::endl;   
      }

      fs<<std::endl;
      fs.close();
    }

    std::string createProcessFile(const char* fname, int irank) {
      std::string file(fname);       
      char temp[FILE_NAME_SIZE];
      snprintf(temp, sizeof(temp), "_%d.dat", irank);  
      file.append(temp);
      return file;
    }
    
    template<typename Vec> 
    inline void printCoords(Vec const& v1, Vec const& v2, Vec const& v3, size_t size, const char* file, int irank) {
      std::fstream fs;
      fs.open(createProcessFile(file, irank).c_str(), std::fstream::out);
      for (size_t i = 0; i < size; ++i) {
        fs << std::setprecision(alogger::precis) << v1[i] << " " << v2[i] << " " << v3[i] << std::endl;
      }
      fs<<std::endl;
      fs.close();
    }

  }
}

#endif
