#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <iostream>
#include <stdint.h>
#include <fstream>
#include "fbind.h"
#include "utils.h"
#include "logger.h"
#include "scattered_field.h"
#include "field_discretizer.h"
#include "source_field.h"

using namespace bemfmm;

int main(int argc, char** argv) {
  int mpisize, mpirank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  bemfmm::node_id = mpirank;  
  alogger::verbose = mpirank==0;
  if(mpirank == 0){
    std::cout << "Running experiment on: " << mpisize << std::endl;
    std::cout << "Underlying args: " << std::endl;
    for(int i = 1; i < argc; ++i) { 
      std::cout << argv[i] << " ";
    }
    std::cout << std::endl;
  }
  alogger::startTimer("Initializing Input");    
  common_data & comm_data = common_data::initCommonData(argc, argv, mpirank, mpisize);
  integral_data const&  int_data = integral_data::initIntegrationData(comm_data);
  alogger::stopTimer("Initializing Input");      
  comm_data.printParams(alogger::stringLength, alogger::verbose); 
  alogger::printTitle("C++ - Wave Scattering");  
  d_complex_t_vec zzsparse;
  d_complex_t_vec diagonal;
  alogger::startTimer("Total Time");  
#if !USE_FMM 
  if(mpisize > 1) {
    if(mpirank == 0) 
      std::cerr << "Dense version is not supported in distributed memory" << std::endl;
    exit(0);
  }
  alogger::startTimer("Building Matrix");
  field_discretizer::discretize(int_data, comm_data, zzsparse, diagonal);  
  alogger::stopTimer("Building Matrix");
#endif      
  alogger::printTime("Initializing Input");
  alogger::startTimer("Fetching Coordinates");
#if USE_PART
  size_t const n = comm_data.nipp * comm_data.init_part_size;  
  d_vector x(n); d_vector y(n); d_vector z(n);
  utils::fetchCoords(comm_data.ntriangle, comm_data.nipp, comm_data.nnpt, 
                 int_data.alphas, int_data.betas, int_data.gammas, comm_data, x, y, z);   
#else
  size_t const n = comm_data.nipp * comm_data.ntriangle;  
  d_vector x(n); d_vector y(n); d_vector z(n);
  utils::fetchCoords(comm_data.ntriangle, comm_data.nipp, comm_data.nnpt, comm_data.nsupan, 
                 int_data.alphas, int_data.betas, int_data.gammas, comm_data.sunod, x, y, z); 
#endif
  //if(comm_data.writeTimingOutputs) utils::printCoords(x, y, z, n, mpirank); 
  alogger::stopTimer("Fetching Coordinates");  
  alogger::startTimer("Solving AX=B");
  int32_vec out_patches; int16_vec out_pt_locs; int out_size;
  source_field source_f(int_data, comm_data, x, y, z);
  d_complex_t_vec rj = source_f.computeSourceField(zzsparse, out_size, out_patches, out_pt_locs);
  alogger::stopTimer("Solving AX=B");    
#if FIND_SCATTERED_FIELD
  scattered_field field_finder(int_data, comm_data,rj, comm_data.nsupan);
  alogger::startTimer("Far Scattered Field");
  field_finder.calculateScatteredField(FAR, mpirank, mpisize, out_size, comm_data.writeTimingOutputs, out_patches, out_pt_locs, "far_sca_soft");
  alogger::stopTimer("Far Scattered Field");
  alogger::startTimer("Near Scattered Field - Far Scheme");
  field_finder.calculateScatteredField(NEAR_FAR_SCHEME, mpirank, mpisize, out_size, comm_data.writeTimingOutputs, out_patches, out_pt_locs, "near_sca_soft_far_scheme");
  alogger::stopTimer("Near Scattered Field - Far Scheme");
  alogger::startTimer("Near Scatterd Field - Near Scheme");
  field_finder.calculateScatteredField(NEAR_NEAR_SCHEME, mpirank, mpisize, out_size, comm_data.writeTimingOutputs, out_patches, out_pt_locs, "far_sca_soft_far_scheme");
  alogger::stopTimer("Near Scatterd Field - Near Scheme");
#endif
  alogger::stopTimer("Total Time");
  alogger::resetTimer();
  MPI_Finalize();
  return 0;
}
