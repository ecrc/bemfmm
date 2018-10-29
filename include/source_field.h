#ifndef SOURCE_FIELD_H
#define SOURCE_FIELD_H
#include "fbind.h"
#include "utils.h"
#include "gmres.h"
#include "fmm/acoustics_wrapper.h"
#include <fstream>

namespace bemfmm {
  class source_field
  {
  private:
    integral_data const& int_data;
    common_data & comm_data;     
    int** const& nsupan;
  	int const& ntriangle;
  	int const& nnpt;
  	int const& nipp;
  	int_data_type const& alphao;
  	int_data_type const& betao;
  	int_data_type const& gammao;
    d_complex_t const& cjvk;
    d_vector& xb;
    d_vector& yb;
    d_vector& zb;   
    double vkk[3];
    double vhh[3];
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
    source_field(integral_data const& _int_data, common_data & _comm_data, 
                d_vector& x, d_vector& y, d_vector& z): 
                int_data(_int_data), comm_data(_comm_data), nsupan(comm_data.nsupan), 
                ntriangle(comm_data.ntriangle),
                nnpt(comm_data.nnpt), nipp(comm_data.nipp),
                alphao(int_data.alphao), betao(int_data.betao), gammao(int_data.gammao),
                cjvk(comm_data.cjvk), xb(x), yb(y), zb(z) {
      for (int i = 0; i < 3; ++i) {
        vkk[i] = comm_data.vkk[i];
        vhh[i] = comm_data.vhh[i];
      }
    }

#if USE_PART
    d_complex_t_vec computeSourceField(d_complex_t_vec const& zzsparse, int& out_size, int32_vec& patches, int16_vec& pointlocs) {
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
      d_complex_t_vec rj(out_size,0.0);
      d_complex_t_vec wb(out_size);     
      for (int i = 0; i < out_size; ++i) {
        wb[i] = ws[pointlocs[i]] * 0.5 / (4.0 * M_PI);
      }
      my_gmres(out_size, ntriangle, nipp, comm_data.nitermax, comm_data.precis, crhs, rj, wb, meta_data, comm_data.fmmVerbose, patches, pointlocs, zzsparse, comm_data.mpirank);      
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
    d_complex_t_vec computeSourceField(d_complex_t_vec const& zzsparse, int& out_size, int32_vec& patches, int16_vec& pointlocs) {
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
      d_complex_t_vec rj(out_size,0.0);
      d_complex_t_vec wb(out_size);     
      for (int i = 0; i < out_size; ++i) {
        wb[i] = ws[pointlocs[i]] * 0.5 / (4.0 * M_PI);
      }
      my_gmres(out_size, ntriangle, nipp, comm_data.nitermax, comm_data.precis, crhs, rj, wb, meta_data, comm_data.fmmVerbose, patches, pointlocs, zzsparse, comm_data.mpirank);      
      if(comm_data.writeTimingOutputs) {
        printVecToFile(rj, rj.size(), "", createProcessFile("iterative_solution_rj.dat",comm_data.mpirank));
      }
      return rj;
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
      return crhs;
#endif
    }
#endif
   };
 }
#endif
