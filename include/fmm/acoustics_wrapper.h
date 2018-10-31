/**
 *
 * @file acoustics_wrapper.h
 *
 * @copyright 2018 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 *
 * @author Mustafa Abduljabbar [mustafa.abduljabbar@kaust.edu.sa] and Mohammed Al Farhan [mohammed.farhan@kaust.edu.sa].
 *
 **/

#ifndef ACOUSTICS_WRAPPER
#define ACOUSTICS_WRAPPER

#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include <fstream>
namespace exafmm {
  struct fmm_data {
    public:
      double kreal;
      double kimag;
      int ncrit; 
      int nthreads; 
      int listbased;
      int tspawn;
      std::string partitioning;
  };
  vec3 cycles;
  Bodies buffer;
  Bounds globalBounds;
  Bodies bbodies;
  Bodies vbodies;
  Bodies fbbodies;
  Bodies fvbodies;
  Cells bcells, fbcells;
  Cells vcells, fvcells;
  int ntriangles;

  BaseMPI * baseMPI;
  BoundBox * boundBox;
  BuildTree * localTree, * fineLocalTree, * globalTree;
  Args* args;
  Partition * partition;
  Traversal * traversal;
  TreeMPI * treeMPI;
  UpDownPass * upDownPass;
  UpDownPass * upDownPass2;
  std::vector<std::vector<double> > nearGauss;
  std::vector<int> patches;
  bool init;
  bool init_low;
  void log_initialize() {
    args->verbose &= baseMPI->mpirank == 0;
    logger::verbose = args->verbose;
    logger::printTitle("FMM Parameters");
    args->print(logger::stringLength, P);
    logger::printTitle("FMM Profiling");
    logger::resetTimer();
    logger::startTimer("Total FMM");
    logger::startPAPI();
  }

  void log_finalize() {
    logger::stopPAPI();
    logger::stopTimer("Total FMM");
    logger::printTitle("Total runtime");
    logger::printTime("Total FMM");
  }
#if USE_PART
    /**
   * initializes FMM
   * @param eps2 is the FMM epsilon 
   * @param fmmAttributes FMM performance control attributes
   * @param nb is the intitial number of bodies
   * @param xb a vector of X coordinates of bodies
   * @param yb a vector of Y coordinates of bodies
   * @param zb a vector of Z coordinates of bodies
   * @param xt a vector of X coordinates of triangular points
   * @param yt a vector of Y coordinates of triangular points
   * @param zt a vector of Z coordinates of triangular points
   * @param patchids a vector of Mesh node ids corresponding to each particle (to avoid self interaction among nodes)
   * @param nearGaussPoints Guass quadrature elements for near-field treatments (when enabled)
   * @param nhdgqp count of Gauss refinements for near-field treatments (when enabled)
   * @param ntriangles_ number of triangular patches
   * @param nipp_ number of integration points per patch
   * @param nearpd near patch distance threshold (to apply p-refinement)
   * @param ws integration points 
   * @param ipolator_near basis vector
  */  
  void FMM_Init(double eps2, fmm_data fmmAttributes, 
           int nb, std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, 
           std::vector<double>& xt, std::vector<double>& yt, std::vector<double>& zt, std::vector<int>& patchids, 
           std::vector<std::vector<double> > nearGaussPoints, int nhdgqp, int ntriangles_,
           int nipp_, double nearpd, std::vector<double> ws, 
           std::vector<std::vector<double> > ipolator_near) {
    int nspawn = fmmAttributes.tspawn;
    const int images = 0;
    const double theta = 0.4;
    const bool useRmax = false;
    const bool useRopt = false;
    const bool verbose = false;  
    init = true;
    init_low = true;
    num_threads(fmmAttributes.nthreads);
    kernel::eps2 = eps2;
    kernel::wavek = complex_t(fmmAttributes.kreal, fmmAttributes.kimag);
    kernel::nhdgqp = nhdgqp;
    kernel::nipp = nipp_;
    ntriangles = ntriangles_;
  #if EXAFMM_NEARF_TREATMENT
    kernel::nearpd = nearpd;
    nearGauss = nearGaussPoints;
  #if EXAFMM_SINGLE
    kernel::ws.resize(ws.size());
    std::copy(ws.begin(), ws.end(), kernel::ws.begin());
    kernel::ipolator_near = ipolator_near;
  #else
    kernel::ws = ws;
    kernel::ipolator_near = ipolator_near;
  #endif
  #endif
    kernel::setup();
    args = new Args;
    baseMPI = new BaseMPI;
    boundBox = new BoundBox(nspawn);
    localTree = new BuildTree(fmmAttributes.ncrit, nspawn);
    fineLocalTree = new BuildTree(64, 1024);
    globalTree = new BuildTree(1, nspawn);
    partition = new Partition(baseMPI->mpirank, baseMPI->mpisize);
    traversal = new Traversal(nspawn, images);
    treeMPI = new TreeMPI(baseMPI->mpirank, baseMPI->mpisize, images);
    upDownPass = new UpDownPass(theta, useRmax, useRopt);
    upDownPass2 = new UpDownPass(1, useRmax, useRopt);
    args->ncrit = fmmAttributes.ncrit;
    args->write = 1;
    args->threads = fmmAttributes.nthreads;
    args->distribution = "external";
    args->dual = (fmmAttributes.listbased == 0);
    args->graft = 0;
    args->images = images;
    args->mutual = 0;
    args->numBodies = 0;
    args->useRopt = useRopt;
    args->nspawn = nspawn;
    args->theta = theta;
    args->partitioning = fmmAttributes.partitioning.c_str();
    //args->verbose = verbose;// & (baseMPI->mpirank == 0);
    args->verbose = verbose & (baseMPI->mpirank == 0);
    //args->verbose = 1;
    args->useRmax = useRmax;
    logger::verbose = args->verbose;
    if(nb <= 0) return; 
    assert((nb % kernel::nipp) == 0);
    vbodies.resize(nb);
    bbodies.resize(nb);
    patches.resize(nb);
    int offset = (ntriangles/baseMPI->mpisize)*baseMPI->mpirank;
  #pragma omp parallel for
    for(int i = 0; i < nb; ++i){
      B_iter B = bbodies.begin() + i;
      int point_index = i;
      B->X[0] = xb[point_index];
      B->X[1] = yb[point_index];
      B->X[2] = zb[point_index];
      
      B->TRI_POINT[0] = xt[point_index];
      B->TRI_POINT[1] = yt[point_index];
      B->TRI_POINT[2] = zt[point_index];

      int patch = patchids[point_index];
      patches[i]=patch;
      assert(patch < ntriangles);
      B->PATCH = patch;
      B->POINT_LOC = point_index%kernel::nipp;
      int relative_patch = patch - offset;
  #if EXAFMM_NEARF_TREATMENT    
      for (int j = 0; j < nhdgqp; ++j) { 
        B->GAUSS_NEAR[j][0] = nearGaussPoints[relative_patch*nhdgqp+j][0];
        B->GAUSS_NEAR[j][1] = nearGaussPoints[relative_patch*nhdgqp+j][1];
        B->GAUSS_NEAR[j][2] = nearGaussPoints[relative_patch*nhdgqp+j][2]; 
      }
  #endif
      B->WEIGHT = 1;
    }
  #pragma omp parallel for
    for(int i = 0; i < nb; ++i){
      B_iter B = vbodies.begin() + i;
      int point_index = i;
      B->X[0] = xb[point_index];
      B->X[1] = yb[point_index];
      B->X[2] = zb[point_index];
      
      B->TRI_POINT[0] = xt[point_index];
      B->TRI_POINT[1] = yt[point_index];
      B->TRI_POINT[2] = zt[point_index];

      int patch = patchids[point_index];
      B->PATCH = patch;
      B->POINT_LOC = point_index%kernel::nipp;
      int relative_patch = patch - offset;
  #if EXAFMM_NEARF_TREATMENT    
      for (int j = 0; j < nhdgqp; ++j) { 
        B->GAUSS_NEAR[j][0] = nearGaussPoints[relative_patch*nhdgqp+j][0];
        B->GAUSS_NEAR[j][1] = nearGaussPoints[relative_patch*nhdgqp+j][1];
        B->GAUSS_NEAR[j][2] = nearGaussPoints[relative_patch*nhdgqp+j][2]; 
      }
  #endif
      B->WEIGHT = 1;
    }
    log_initialize();  
  }

 /**
   * partitions FMM
   * @param nb is the intitial number of bodies (input/output)
   * @param xb a vector of X coordinates of bodies (input/output)
   * @param yb a vector of Y coordinates of bodies (input/output)
   * @param zb a vector of Z coordinates of bodies (input/output)
   * @param xt a vector of X coordinates of triangular nodes (input/output)
   * @param yt a vector of Y coordinates of triangular nodes (input/output)
   * @param zt a vector of Z coordinates of triangular nodes (input/output)
   * @param patch a vector of Mesh node ids corresponding to each particle to avoid self interaction among nodes (input/output)
   * @param loc the location of eatch particle with the triangular mesh (input/output)
  */  
  void FMM_Partition(int& nb, std::vector<double>&xb, std::vector<double>&yb, std::vector<double>&zb, std::vector<double>&xt, std::vector<double>&yt, std::vector<double>&zt, std::vector<int>& patch, std::vector<short>& loc) { 
   logger::printTitle("Partition Profiling");
   Bounds localBounds;
   if(nb > 0) {
      localBounds = boundBox->getBounds(bbodies);
      localBounds = boundBox->getBounds(vbodies, localBounds);
    } else {
      localBounds.Xmin = 0;
      localBounds.Xmax = 0;
    }
    
    globalBounds = baseMPI->allreduceBounds(localBounds);
    cycles = globalBounds.Xmax - globalBounds.Xmin;
    if(baseMPI->mpisize == 1)  { 
      for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) { 
        int i = B-bbodies.begin();
        B->IBODY = i;
      }
      for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
        int i = B-vbodies.begin();
        B->IBODY = i;
      }
      nb = bbodies.size();
      return;
    }
    int nPatches = bbodies.size()/kernel::nipp;
    assert(bbodies.size()%kernel::nipp == 0);
    Bodies partitionBodies(nPatches);
    for (B_iter B=bbodies.begin(), BP=partitionBodies.begin(); B!=bbodies.end(); B+=kernel::nipp, BP++) {
       BP->X[0] = B->X[0];
       BP->X[1] = B->X[1];
       BP->X[2] = B->X[2];
       BP->PATCH = B->PATCH;
       BP->WEIGHT = 1;
       BP->IBODY = B - bbodies.begin();
    }
    partition->partition(partitionBodies, globalBounds, args->partitioning);
    Bodies tempBodies(bbodies);
    int currPatch = 0;
    B_iter bbodiesBegin = bbodies.begin();
    B_iter vbodiesBegin = vbodies.begin();
    for (B_iter B=partitionBodies.begin(); B!=partitionBodies.end(); B++) { 
      B_iter TB = tempBodies.begin() + B->IBODY;
      B_iter BB = bbodiesBegin+currPatch*kernel::nipp;
      B_iter VB = vbodiesBegin+currPatch*kernel::nipp;
      for (int i = 0; i < kernel::nipp; ++i) {
        assert(B->IRANK >= 0 && B->IRANK < baseMPI->mpisize);
        (TB+i)->IRANK = B->IRANK;
        *(BB+i) = *(TB+i);
        *(VB+i) = *(TB+i);
      }    
      currPatch++;
    }
    assert(bbodiesBegin+currPatch*kernel::nipp == bbodies.end());
    nb = bbodies.size();
    int global; 
    MPI_Allreduce(&nb, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(global != (ntriangles * kernel::nipp)) std::cout << global  << "  : " << ntriangles*kernel::nipp << std::endl;
    assert(global == (ntriangles * kernel::nipp));           
    bbodies = treeMPI->commBodies(bbodies);
    vbodies = treeMPI->commBodies(vbodies);  
    nb = bbodies.size();
    MPI_Allreduce(&nb, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(global != (ntriangles * kernel::nipp)) std::cout << global  << "  : " << ntriangles*kernel::nipp << std::endl;
    assert(global == (ntriangles * kernel::nipp));        
    patch.resize(nb);
    loc.resize(nb);
    xb.resize(nb);
    yb.resize(nb);
    zb.resize(nb);
    xt.resize(nb);
    yt.resize(nb);
    zt.resize(nb);
    for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
      int i = B-bbodies.begin();
      patch[i] = B->PATCH;
      assert(patch[i] < ntriangles);
      loc[i] = B->POINT_LOC;
      xb[i] = B->X[0];
      yb[i] = B->X[1];
      zb[i] = B->X[2];
      xt[i] = B->TRI_POINT[0];
      yt[i] = B->TRI_POINT[1];
      zt[i] = B->TRI_POINT[2];
      B->IBODY = i;
    }
    for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
      int i = B-vbodies.begin();
      B->IBODY = i;
    }
  }
#else
  /**
   * initializes FMM
   * @param eps2 is the FMM epsilon 
   * @param fmmAttributes FMM performance control attributes
   * @param nb is the intitial number of bodies
   * @param xb a vector of X coordinates of bodies
   * @param yb a vector of Y coordinates of bodies
   * @param zb a vector of Z coordinates of bodies
   * @param patchids a vector of Mesh node ids corresponding to each particle (to avoid self interaction among nodes)
   * @param nearGaussPoints Guass quadrature elements for near-field treatments (when enabled)
   * @param nhdgqp count of Gauss refinements for near-field treatments (when enabled)
   * @param ntriangles_ number of triangular patches
   * @param nipp_ number of integration points per patch
   * @param nearpd near patch distance threshold (to apply p-refinement)
   * @param ws integration points 
   * @param ipolator_near basis vector
  */  
  void FMM_Init(double eps2, fmm_data fmmAttributes, 
           int nb, std::vector<double>& xb, std::vector<double>& yb, std::vector<double>& zb, std::vector<int>& patchids, 
                           std::vector<std::vector<double> > nearGaussPoints, int nhdgqp, int ntriangles_,
                           int nipp_, double nearpd, std::vector<double> ws, 
                           std::vector<std::vector<double> > ipolator_near) {
    int nspawn = fmmAttributes.tspawn;
    const int images = 0;
    const double theta = 0.4;
    const bool useRmax = false;
    const bool useRopt = false;
    const bool verbose = false;  
    init = true;
    init_low = true;
    num_threads(fmmAttributes.nthreads);
    kernel::eps2 = eps2;
    kernel::wavek = complex_t(fmmAttributes.kreal, fmmAttributes.kimag);
    kernel::nhdgqp = nhdgqp;
    kernel::nipp = nipp_;
    ntriangles = ntriangles_;
  #if EXAFMM_NEARF_TREATMENT
    kernel::nearpd = nearpd;
    nearGauss = nearGaussPoints;
  #if EXAFMM_SINGLE
    kernel::ws.resize(ws.size());
    std::copy(ws.begin(), ws.end(), kernel::ws.begin());
    kernel::ipolator_near = ipolator_near;
  #else
    kernel::ws = ws;
    kernel::ipolator_near = ipolator_near;
  #endif
  #endif
    kernel::setup();
    args = new Args;
    baseMPI = new BaseMPI;
    boundBox = new BoundBox(nspawn);
    localTree = new BuildTree(fmmAttributes.ncrit, nspawn);
    fineLocalTree = new BuildTree(64, 1024);
    globalTree = new BuildTree(1, nspawn);
    partition = new Partition(baseMPI->mpirank, baseMPI->mpisize);
    traversal = new Traversal(nspawn, images);
    treeMPI = new TreeMPI(baseMPI->mpirank, baseMPI->mpisize, images);
    upDownPass = new UpDownPass(theta, useRmax, useRopt);
    upDownPass2 = new UpDownPass(1, useRmax, useRopt);
    args->ncrit = fmmAttributes.ncrit;
    args->write = 1;
    args->threads = fmmAttributes.nthreads;
    args->distribution = "external";
    args->dual = (fmmAttributes.listbased == 0);
    args->graft = 0;
    args->images = images;
    args->mutual = 0;
    args->numBodies = 0;
    args->useRopt = useRopt;
    args->nspawn = nspawn;
    args->theta = theta;
    args->partitioning = fmmAttributes.partitioning.c_str();
    //args->verbose = verbose;// & (baseMPI->mpirank == 0);
    args->verbose = verbose & (baseMPI->mpirank == 0);
    //args->verbose = 1;
    args->useRmax = useRmax;
    logger::verbose = args->verbose;
    if(nb <= 0) return; 
    assert((nb % kernel::nipp) == 0);
    vbodies.resize(nb);
    bbodies.resize(nb);
    patches.resize(nb);
#pragma omp parallel for
    for(int i = 0; i < nb; ++i){
      B_iter B = bbodies.begin() + i;
      int point_index = i;
      B->X[0] = xb[point_index];
      B->X[1] = yb[point_index];
      B->X[2] = zb[point_index];

      int patch = patchids[point_index];
      patches[i]=patch;
      B->PATCH = patch;
      B->POINT_LOC = point_index%kernel::nipp;
  #if EXAFMM_NEARF_TREATMENT    
      for (int j = 0; j < nhdgqp; ++j) { 
        B->GAUSS_NEAR[j][0] = nearGaussPoints[patch*nhdgqp+j][0];
        B->GAUSS_NEAR[j][1] = nearGaussPoints[patch*nhdgqp+j][1];
        B->GAUSS_NEAR[j][2] = nearGaussPoints[patch*nhdgqp+j][2]; 
      }
  #endif
      B->WEIGHT = 1;
    }
#pragma omp parallel for 
    for(int i = 0; i < nb; ++i){
      B_iter B = vbodies.begin() + i;
      int point_index = i;
      B->X[0] = xb[point_index];
      B->X[1] = yb[point_index];
      B->X[2] = zb[point_index];

      int patch = patchids[point_index];
      B->PATCH = patch;
      B->POINT_LOC = point_index%kernel::nipp;
  #if EXAFMM_NEARF_TREATMENT    
      for (int j = 0; j < nhdgqp; ++j) { 
        B->GAUSS_NEAR[j][0] = nearGaussPoints[patch*nhdgqp+j][0];
        B->GAUSS_NEAR[j][1] = nearGaussPoints[patch*nhdgqp+j][1];
        B->GAUSS_NEAR[j][2] = nearGaussPoints[patch*nhdgqp+j][2]; 
      }
  #endif
      B->WEIGHT = 1;
    }
    log_initialize();  
  }
  /**
   * partitions FMM
   * @param nb is the intitial number of bodies (input/outpu)
   * @param xb a vector of X coordinates of bodies (input/outpu)
   * @param yb a vector of Y coordinates of bodies (input/outpu)
   * @param zb a vector of Z coordinates of bodies (input/outpu)
   * @param patch a vector of Mesh node ids corresponding to each particle to avoid self interaction among nodes (input/output)
   * @param loc the location of eatch particle with the triangular mesh (input/output)
  */  
  void FMM_Partition(int& nb, std::vector<double>&xb, std::vector<double>&yb, std::vector<double>&zb, std::vector<int>& patch, std::vector<short>& loc) { 
   logger::printTitle("Partition Profiling");
   Bounds localBounds;
   if(nb > 0) {
      localBounds = boundBox->getBounds(bbodies);
      localBounds = boundBox->getBounds(vbodies, localBounds);
    } else {
      localBounds.Xmin = 0;
      localBounds.Xmax = 0;
    }
    
    globalBounds = baseMPI->allreduceBounds(localBounds);
    cycles = globalBounds.Xmax - globalBounds.Xmin;
    if(baseMPI->mpisize == 1)  { 
      for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) { 
        int i = B-bbodies.begin();
        B->IBODY = i;
      }
      for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
        int i = B-vbodies.begin();
        B->IBODY = i;
      }
      nb = bbodies.size();
      return;
    }
    int nPatches = bbodies.size()/kernel::nipp;
    Bodies partitionBodies(nPatches);
    for (B_iter B=bbodies.begin(), BP=partitionBodies.begin(); B!=bbodies.end(); B+=kernel::nipp, BP++) {
       BP->X[0] = B->X[0];
       BP->X[1] = B->X[1];
       BP->X[2] = B->X[2];
       BP->PATCH = B->PATCH;
       BP->WEIGHT = 1;
       BP->IBODY = B - bbodies.begin();
    }
    partition->partition(partitionBodies, globalBounds,args->partitioning);
    Bodies tempBodies(bbodies);
    int currPatch = 0;
    B_iter bbodiesBegin = bbodies.begin();
    B_iter vbodiesBegin = vbodies.begin();
    for (B_iter B=partitionBodies.begin(); B!=partitionBodies.end(); B++) { 
      B_iter TB = tempBodies.begin() + B->IBODY;
      B_iter BB = bbodiesBegin+currPatch*kernel::nipp;
      B_iter VB = vbodiesBegin+currPatch*kernel::nipp;
      for (int i = 0; i < kernel::nipp; ++i) {
        (TB+i)->IRANK = B->IRANK;
        *(BB+i) = *(TB+i);
        *(VB+i) = *(TB+i);
      }    
      currPatch++;
    }
    bbodies = treeMPI->commBodies(bbodies);
    vbodies = treeMPI->commBodies(vbodies);
    nb = bbodies.size();
    patch.resize(nb);
    loc.resize(nb);
    xb.resize(nb);
    yb.resize(nb);
    zb.resize(nb);
    for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
      int i = B-bbodies.begin();
      patch[i] = B->PATCH;
      loc[i] = B->POINT_LOC;
      xb[i] = B->X[0];
      yb[i] = B->X[1];
      zb[i] = B->X[2];
      B->IBODY = i;
    }
    for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
      int i = B-vbodies.begin();
      B->IBODY = i;
    }
  }
#endif
  /**
  * Builds FMM Tree
  */  
  void FMM_BuildTree() {
    Bounds localBoundsB = boundBox->getBounds(bbodies);
    bcells = localTree->buildTree(bbodies, buffer, localBoundsB); 
    Bounds localBoundsV = boundBox->getBounds(vbodies);
    vcells = localTree->buildTree(vbodies, buffer, localBoundsV);
  }
  /**
   * Applies FMM
   * @param vi target values (output)
   * @param vb source values (input)
   * @param wb source weights (input)
   * @param verbose output verbosity (input)
  */  
  void FMM_MatVec(std::vector<std::complex<double> >& vi, std::vector<std::complex<double> >const& vb, std::vector<std::complex<double> >const& wb, bool verbose) { 
    args->verbose = verbose;
    log_initialize();  
    for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
      B->SRC    = 1.0;
      B->QWEIGHT = 1.0;
      B->TRG    = 0;
      B->ICELL = 0;
    }
    for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
      B->SRC    = vb[B->IBODY];
      B->QWEIGHT = wb[B->IBODY];
      B->TRG    = 0;
      B->ICELL = 0;
    }    
    if(init) {
      FMM_BuildTree();
      upDownPass->upwardPass(bcells);
      upDownPass->upwardPass(vcells);
      init = 0;
    } else { 
      upDownPass2->upwardPass(bcells, false);
      upDownPass2->upwardPass(vcells, false);
    }        
      treeMPI->setLET(vcells, cycles);
  #pragma omp parallel sections
  { 
  #pragma omp section
  {
  #if EXAFMM_USE_DISTGRAPH 
      treeMPI->initDistGraph(globalBounds);       
      treeMPI->commDistGraph(cycles);   
  #else
      treeMPI->commBodies();
      treeMPI->commCells();
  #endif       
  } 
  #pragma omp section
  {
      traversal->initListCount(bcells);
      traversal->initWeight(bcells);
      traversal->traverse(bcells, vcells, cycles, args->dual, args->mutual);
  }
  }
    if (baseMPI->mpisize > 1) {
      if (args->graft) {
        treeMPI->linkLET();
        Bodies gbodies = treeMPI->root2body();
        Cells jcells = globalTree->buildTree(gbodies, buffer, globalBounds);
        treeMPI->attachRoot(jcells);
        traversal->traverse(bcells, jcells, cycles, args->dual, false);
      } else {
        for (int irank=0; irank<baseMPI->mpisize; irank++) {
    Cells jcells;
    treeMPI->getLET(jcells, (baseMPI->mpirank+irank)%baseMPI->mpisize);
    traversal->traverse(bcells, jcells, cycles, args->dual, false);
        }
      }
    }

    upDownPass->downwardPass(bcells);
    if(verbose) {
      localTree->printTreeData(bcells); 
      traversal->printTraversalData();
      traversal->writeTraversalData(baseMPI->mpirank, bbodies.size());
      logger::writeTime(baseMPI->mpirank);
    }
    log_finalize();
    for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
      assert(B->IBODY < bbodies.size() && B->IBODY >= 0);
      vi[B->IBODY] = B->TRG[0];
    }
  }
  /**
   * Finalizes FMM
  */  
  void FMM_Finalize() {
    delete args;
    delete baseMPI;
    delete boundBox;
    delete localTree;
    delete fineLocalTree;
    delete globalTree;
    delete partition;
    delete traversal;
    delete treeMPI;
    delete upDownPass;
  }
  /**
   * calculates result using N-Body direct fomr the first nb bodies
   * @param nb number of bodies
   * @param vi target values
   * @param addresses locations of sample bodies
  */  
  void Direct(int nb, std::complex<double>* vi, std::vector<int>& addresses) {
    Bodies ibodies(bbodies.begin(), bbodies.begin() + nb);
    Bodies jbodies(vbodies.begin(), vbodies.end());
    for (B_iter B=ibodies.begin(); B!=ibodies.end(); B++) {      
      B->TRG    = 0;
      B->ICELL = 0;
    }
    for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {      
      B->TRG    = 0;
      B->ICELL = 0;
    }

    for (int irank=0; irank<baseMPI->mpisize; irank++) {
      if (args->verbose) std::cout << "Direct loop          : " << irank+1 << "/" << baseMPI->mpisize << std::endl;
      treeMPI->shiftBodies(jbodies);
      traversal->direct(ibodies, jbodies, cycles);
    }
    for (B_iter B=ibodies.begin(); B!=ibodies.end(); B++) {
      int i = B-ibodies.begin();
      addresses[i] = B->IBODY;
      vi[B->IBODY] = B->TRG[0];
    }
  }  

  /**
   * calculates result using N-Body direct fomr the sample nb bodies
   * @param nb number of sample bodies
   * @param vi target values
   * @param addresses locations of sample bodies
  */  
  void DirectSample(int nb, std::vector<std::complex<double> >& vi, std::vector<int>& addresses) {
    Bodies ibodies(nb);
    Bodies jbodies(vbodies.begin(), vbodies.end());
    srand(time(NULL));
    for (B_iter B=ibodies.begin(); B!=ibodies.end(); B++) {      
      *B = bbodies[rand()%bbodies.size()];
      B->TRG    = 0;
      B->ICELL = 0;
    }
    for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {      
      B->TRG    = 0;
      B->ICELL = 0;
    }

    for (int irank=0; irank<baseMPI->mpisize; irank++) {
      if (args->verbose) std::cout << "Direct loop          : " << irank+1 << "/" << baseMPI->mpisize << std::endl;
      treeMPI->shiftBodies(jbodies);
      traversal->direct(ibodies, jbodies, cycles);
    }
    for (B_iter B=ibodies.begin(); B!=ibodies.end(); B++) {
      int i = B-ibodies.begin();
      addresses[i] = B->IBODY;
      vi[B->IBODY] = B->TRG[0];
    }
  }

}

#endif
