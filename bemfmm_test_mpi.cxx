#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <iostream>
#include <stdint.h>
#include <fstream>
#include "fbind.h"
#include "utils.h"
#include "BEMSolver.h"

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
  BEMSolver solver(int_data, comm_data);
  alogger::stopTimer("Initializing Input");      
  comm_data.printParams(alogger::stringLength, alogger::verbose); 
  alogger::printTitle("C++ - Wave Scattering");  
  alogger::startTimer("Total Time");  
#if !USE_FMM 
  if(mpisize > 1) {
    if(mpirank == 0) 
      std::cerr << "Dense version is not supported in distributed memory" << std::endl;
    exit(0);
  }
  alogger::startTimer("Building Matrix");
  solver.buildDenseMatrix();
  alogger::stopTimer("Building Matrix");
#endif      
  alogger::printTime("Initializing Input"); 
  alogger::startTimer("Solving AX=B");
  int32_vec out_patches; int16_vec out_pt_locs; int out_size;
  d_complex_t_vec rj = solver.computeSourceField(out_size, out_patches, out_pt_locs);
  alogger::stopTimer("Solving AX=B");    
  alogger::startTimer("Far Scattered Field");
  solver.calculateScatteredField(FAR, out_size, comm_data.writeTimingOutputs, out_patches, out_pt_locs, "scatterd_field");
  alogger::stopTimer("Far Scattered Field");
  alogger::startTimer("Near Scattered Field - Far Scheme");
  solver.calculateScatteredField(NEAR_FAR_SCHEME, out_size, comm_data.writeTimingOutputs, out_patches, out_pt_locs, "near_field_scheme1");
  alogger::stopTimer("Near Scattered Field - Far Scheme");
  alogger::startTimer("Near Scatterd Field - Near Scheme");
  solver.calculateScatteredField(NEAR_NEAR_SCHEME, out_size, comm_data.writeTimingOutputs, out_patches, out_pt_locs, "near_field_scheme2");
  alogger::stopTimer("Near Scatterd Field - Near Scheme");
  alogger::stopTimer("Total Time");
  alogger::resetTimer();
  MPI_Finalize();
  return 0;
}