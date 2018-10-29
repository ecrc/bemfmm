#ifndef PARMETIS_WRAPPER_H
#define PARMETIS_WRAPPER_H

namespace parmetis_wrapper{
#include "parmetis.h"
#include <iostream>
  void partitionCoodsParmetis(size_t* counts,float* xyz, int* part, int* ndims, int p,  MPI_Comm* comm) {
    idx_t* vtxdist = new idx_t[p+1];
    vtxdist[0] = 0;
    for(int i = 1; i <= p; ++i){
      vtxdist[i] = vtxdist[i-1] + counts[i-1];
    }
    if(ParMETIS_V3_PartGeom(vtxdist, ndims , xyz, part, comm) == METIS_ERROR){
      std::cout << "Error in METIS Call" << std::endl;
    } 
    delete[] vtxdist;
  }
}

#endif
