/**
 *
 * @file global_data.h
 *
 * @copyright 2018 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 *
 * @author Mustafa Abduljabbar [mustafa.abduljabbar@kaust.edu.sa] and Mohammed Al Farhan [mohammed.farhan@kaust.edu.sa].
 *
 **/

#ifndef GLOBAL_DATA
#define GLOBAL_DATA

#include <assert.h>
#include <fstream>
#include <sstream>
#include "logger.h"
#include "utils.h"
#include "acoustics.h"

namespace bemfmm {
  //!  common_data class
  /*!
    Holds the data needed throughout the execution of the solver
  */     
  class common_data {
 public:
    //! A public variable.
    /*!
      Wave constant
    */
    d_complex_t const cjj;
    //! A public variable.
    /*!
      Number of Integral Points per Patch.
    */
    int nipp;
    //! A public variable.
    /*!
      Number of High degree Gauss Quadrature Points.
    */
    int nhdgqp;
    //! A public variable.
    /*!
      Number of Mfie line integral points
    */
    int nlqp;
     //! A public variable.
    /*!
      Number of Gauss points for RCS field
    */
    int ngprcs;
    //! A public variable.
    /*!
      Number of efie contour integral points
    */
    int nclgp;
    //! A public variable.
    /*!
      0: read ANSYS file;  1: read IDEA-s file
    */
    int geomFormat;
    //! A public variable.
    /*!
      Number of triangular patches
    */
    int ntriangle;
    //! A public variable.
    /*!
      Number of nodes per triangular patch
    */
    int nnpt;
    //! A public variable.
    /*!
      Number of RCS phi's
    */
    int nphi_rcs;
    //! A public variable.
    /*!
      Number of RCS theta's
    */
    int nthe_rcs;    
    //! A public variable.
    /*!
      The starting RCS theta
    */
    double the_sta_degree;
    //! A public variable.
    /*!
      The ending RCS theta
    */    
    double the_end_degree;
    //! A public variable.
    /*!
      The starting RCS phi
    */        
    double phi_sta_degree;
    //! A public variable.
    /*!
      The ending RCS phi
    */            
    double phi_end_degree;
    //! A public variable.
    /*!
      Number of near field observation phi's
    */    
    int nphi_nf;
    //! A public variable.
    /*!
      Number of near field observation phi's
    */        
    int nthe_nf;
    //! A public variable.
    /*!
      Single frequency
    */
    double frequency;
    //! A public variable.
    /*!
      The starting near-field observation theta
    */       
    double the_sta_degree_nf;
    //! A public variable.
    /*!
      The ending near-field observation theta
    */           
    double the_end_degree_nf;
    //! A public variable.
    /*!
      The starting near-field observation phi
    */           
    double phi_sta_degree_nf;
    //! A public variable.
    /*!
      The ending near-field observation phi
    */           
    double phi_end_degree_nf;
    //! A public variable.
    /*!
      Radius of near field observation points (m)
    */   
    double rnfobp;
    //! A public variable.
    /*!
      Wave parameter (k)
    */       
    d_complex_t cjvk;
    //! A public variable.
    /*!
      Number of nodes (DoF)
    */       
    int nnod;
    //! A public variable.
    /*!
       Max distance (m) to apply for P-Mesh refinement 
    */
    double nearpd;
    //! A public variable.
    /*!
      basis vector
    */
    double vkk[3];
    //! A public variable.
    /*!
      basis vector
    */    
    double vee[3];
    //! A public variable.
    /*!
      basis vector
    */
    double vhh[3];

    double the, phi, ethetad, ephid;   
    //! A public variable.
    /*!
       Max number of iterations
    */
    int nitermax;
    //! A public variable.
    /*!
       Precision of iterative solver.
    */
    double precis;
    //! A public variable.
    /*!
      non-zero for detailed FMM timing
    */
    int fmmVerbose;
    //! A public variable.
    /*!
      FMM attributes struct
    */
    exafmm::fmm_data fmmAttributes;
    //! A public variable.
    /*!
      Compare to FMM direct
    */
    bool checkFMMDirect;
    //! A public variable.
    /*!
      GMRES restart iteration
    */
    int gmresRestart;   
    //! A public variable.
    /*!
      Write timing outputs for FMM
    */
    bool writeTimingOutputs;
    //! A public variable.
    /*!
      Curremt MPI rank
    */
    int mpirank;
    //! A public variable.
    /*!
      Current MPI size
    */
    int mpisize;
    //! A public variable.
    /*!
      Geometry input file
    */
    std::string geomatryFileName;    
    //! A public variable.
    /*!
      Coordinates of nodes
    */
    double** sunod;
    //! A public variable.
    /*!
      Indices of nodes
    */
    int** nsupan;
    //! A public variable.
    /*!
      1D array clone of coordinates of nodes
    */    
    double* sunod_;
    //! A public variable.
    /*!
      1D array clone of indices of nodes
    */        
    int* nsupan_;
    //! A public variable.
    /*!
      The Args structure
    */        
    static Args* args;    
    //! A public variable.
    /*!
      The partition size (number of points in this partition)
    */        
    int init_part_size;
#if USE_PART
    //! A public variable.
    /*!
      The start triangle
    */      
    int start_triangle;
    //! A public variable.
    /*!
      The end triangle
    */      
    int end_triangle;    
    //! A public variable.
    /*!
      The x coordinates for local points
    */      
    d_vector xt;
    //! A public variable.
    /*!
      The y coordinates for local points
    */      
    d_vector yt;
    //! A public variable.
    /*!
      The z coordinates for local points
    */          
    d_vector zt;   
#endif
#if USE_PART
    //! A public variable.
    /*!
      The coordinates for triangles
    */         
    double* triangles; 
#endif
    /**
     * initializes the global data for the solver
     * @param argc count passed from main function.
     * @param argv is the arugments string.
     * @param mpirank is the mpi rank.
     * @param mpisize is the mpi size.
    */    
    static common_data& initCommonData(int argc, char** argv, int mpirank, int mpisize) { 
      args = new Args(argc, argv);
      static common_data instance(mpirank, mpisize);
      return instance;
    }
    /**
     * A destructor.
     * destroys geometry data.
    */   
    ~common_data() {
      delete[] sunod_;
      delete[] nsupan_;      
      delete[] sunod;
      delete[] nsupan;
#if USE_PART
      delete triangles;
#endif
      delete common_data::args;
    }
    /**     
     * a print member function 
     * prints solver config params.
    */  
    void printParams(int stringLength, bool verbose) const {
      if (verbose) {
        alogger::printTitle("Wave Scattering Parameters");
        alogger::logItem(stringLength,"ntriangle",ntriangle);
        alogger::logItem(stringLength,"nipp",nipp);
        alogger::logItem(stringLength,"unknowns",ntriangle*nipp);
        alogger::logItem(stringLength,"memory",((pow(ntriangle*nipp,2))*16)/(1<<20));
        alogger::logItem(stringLength,"the_start",the_sta_degree);
        alogger::logItem(stringLength,"the_end",the_end_degree);
        alogger::logItem(stringLength,"phi_start",phi_sta_degree);
        alogger::logItem(stringLength,"phi_end",phi_end_degree);
        alogger::logItem(stringLength,"nphi_rcs",nphi_rcs);
        alogger::logItem(stringLength,"nthe_rcs",nthe_rcs);
        alogger::logItem(stringLength,"the_start_nf",the_sta_degree_nf);
        alogger::logItem(stringLength,"the_end_nf",the_end_degree_nf);
        alogger::logItem(stringLength,"phi_start_nf",phi_sta_degree_nf);
        alogger::logItem(stringLength,"phi_end_nf",phi_end_degree_nf);
        alogger::logItem(stringLength,"nphi_nf",nphi_nf);
        alogger::logItem(stringLength,"nthe_nf",nthe_nf);
        alogger::logItem(stringLength,"nhdgqp",nhdgqp);
        alogger::logItem(stringLength,"nlqp",nlqp);
        alogger::logItem(stringLength,"nclgp",nclgp);
        alogger::logItem(stringLength,"nearpd",nearpd);
        alogger::logItem(stringLength,"ngprcs",ngprcs);
        alogger::logItem(stringLength,"rnfobp",rnfobp);
        alogger::logItem(stringLength,"gmres precis",precis);       
        alogger::logItem(stringLength,"gmres restart",gmresRestart);        
        alogger::logItem(stringLength,"frequency",frequency);
      }
    }

  private:
    /**
     * A private constructor to keep a singleton instance of the class.
     * @param irank is rank.
     * @param isize is size.
    */
    common_data(int irank, int isize) : mpirank(irank), mpisize(isize), cjj(0.0,1.0) {
      const char* inputArgsFile = common_data::args->configfile.c_str();      
      readInputArgs(inputArgsFile);     
      mpirank = irank;
      mpisize = isize;      
      if(geomFormat == 0)
        readTextualAnsysGeom(geomatryFileName.c_str());
      else if(utils::endsWith(geomatryFileName,"bin")){ 
        readBinaryIDEAsGeom(geomatryFileName.c_str());
      }
      else
        readTextualIDEAsGeom(geomatryFileName.c_str());
    }
    /**
     * A private function to broadcast small geometry data.
    */
    void sendRecvGeom(){
      MPI_Bcast(&nnod, 1, MPI_INT,0, MPI_COMM_WORLD);
      MPI_Bcast(&nnpt, 1, MPI_INT,0, MPI_COMM_WORLD);
      MPI_Bcast(&ntriangle, 1, MPI_INT,0, MPI_COMM_WORLD);
      if(mpirank != 0) {
        sunod_ = new double[dim*nnod];
        nsupan_ = new int[ntriangle*nnpt];
        sunod = new double*[dim];
        nsupan = new int*[ntriangle];
      }
      size_t ss = dim*nnod / 4;
      size_t rs = dim*nnod % 4;
      MPI_Bcast(sunod_, ss, MPI_DOUBLE,0, MPI_COMM_WORLD);
      MPI_Bcast(sunod_+ss, ss, MPI_DOUBLE,0, MPI_COMM_WORLD);
      MPI_Bcast(sunod_+2*ss, ss, MPI_DOUBLE,0, MPI_COMM_WORLD);
      MPI_Bcast(sunod_+3*ss, ss, MPI_DOUBLE,0, MPI_COMM_WORLD);
      if(rs > 0){
        MPI_Bcast(sunod_+4*ss, rs, MPI_DOUBLE,0, MPI_COMM_WORLD);
      }
      size_t psize = ntriangle*nnpt;
      if(psize < INT32_MAX){
        MPI_Bcast(nsupan_,ntriangle*nnpt, MPI_INT,0, MPI_COMM_WORLD);        
      } else { 
        int* ns = nsupan_;
        int iters = psize / INT32_MAX;
        int remainder = psize % INT32_MAX;
        for(int i = 0; i < iters; ++i) {
          MPI_Bcast(ns,INT32_MAX, MPI_INT,0, MPI_COMM_WORLD);      
          ns += INT32_MAX;
        }
        if(remainder > 0){ 
          MPI_Bcast(ns,remainder, MPI_INT, 0, MPI_COMM_WORLD);        
        } 
      }
      if(mpirank != 0) {
        for (int i = 0; i < dim; ++i) 
          sunod[i] = (sunod_+i*nnod);        
        for (int i = 0; i < ntriangle; ++i) 
          nsupan[i] = (nsupan_+i*nnpt);          
      }
    }
    /**
     * A private reader for textual geometry data in Ansys format
    */
    void readTextualAnsysGeom(const char* file_name) { 
      if(mpirank == 0) {       
        std::ifstream input_file(file_name);
        std::stringstream ss;
        std::string line;
        int npat, nod;
        if(nnpt == 6) {
          for (int i = 0; i < dim; ++i)
            std::getline(input_file,line);

          input_file >> npat >> nnod;
          std::getline(input_file,line);
          ntriangle=npat;
          init_part_size = ntriangle;
          sunod_ = new double[dim*nnod];
          nsupan_ = new int[npat*nnpt];
          sunod = new double*[dim];
          nsupan = new int*[npat];
          for (int i = 0; i < dim; ++i) 
            sunod[i] = (sunod_+i*nnod);        
          for (int i = 0; i < npat; ++i) 
            nsupan[i] = (nsupan_+i*nnpt);

          int pp = npat % 20;
          int qq = (npat-pp) / 20;  
          int n,tri;

          for (int i = 0; i < qq; ++i) {
            std::getline(input_file,line);
            std::getline(input_file,line);
            std::getline(input_file,line);
            for (int j = 0; j < 20; ++j) {
              std::getline(input_file,line);
              ss.str(line);
              ss >> std::skipws >> tri >> n >> n >> n >> n >> n; 
              for (int k = 0; k < 6; ++k) {
                ss >> nsupan[tri-1][k];
                nsupan[tri-1][k]--;
              }
            }
          }

          if(pp != 0) {
            std::getline(input_file,line);
            std::getline(input_file,line);
            std::getline(input_file,line);
            for (int i = 0; i < pp; ++i) {
              std::getline(input_file,line);
              ss.str(line);
              ss >> std::skipws >> tri >> n >> n >> n >> n >> n; 
              for (int k = 0; k < 6; ++k) {
                ss >> nsupan[tri-1][k];
                nsupan[tri-1][k]--;
              }
            }
          }

          pp = nnod % 20;
          qq = (nnod - pp) / 20;

          for (int i = 0; i < 3; ++i) 
            std::getline(input_file,line);


          for (int i = 0; i < qq; ++i) {
            std::getline(input_file,line);
            std::getline(input_file,line);
            for (int j = 0; j < 20; ++j) {
              std::getline(input_file,line);
              ss.str(line);
              ss >> std::skipws >> nod;
              ss >> sunod[0][nod - 1] >> sunod[1][nod - 1] >> sunod[2][nod - 1];
            }
          }

          if(pp != 0) {
            std::getline(input_file,line);
            std::getline(input_file,line);
            for (int i = 0; i < pp; ++i) {
              std::getline(input_file,line);
              ss.str(line);
              ss >> std::skipws >> nod;
              ss >> sunod[0][nod - 1] >> sunod[1][nod - 1] >> sunod[2][nod - 1];
            }
          }

        } else {
          std::cout<<"cannot read higher order mesh now" << std::endl;
        }
        input_file.close();
      } 
      sendRecvGeom();
    } 
#if USE_PART
      /**
       * A private partitioned reader for binary geometry data in I-DEAs format
      */    
       void readBinaryIDEAsGeom(const char* file_name) { 
         std::ifstream input_file(file_name, std::ios::in|std::ios::binary);
         if(input_file.fail()) { 
           std::cout << "failed to open file " << file_name  <<  " by rank " << mpirank << " make sure file exists and access permissions are ok" <<std::endl;
           exit(EXIT_FAILURE);
         }

         if(nnpt == 6) {  
           char* firstRead = new char[sizeof(int)];
           input_file.seekg (0, std::ios::beg);
           input_file.read(firstRead, sizeof(int));
           input_file.read(firstRead, sizeof(int));
           int npat = *((int*)firstRead);
           ntriangle=npat;
           size_t chunk = ntriangle/mpisize;
           input_file.seekg (sizeof(double)*(chunk*3*6*mpirank+1), std::ios::beg);
           start_triangle = 0;
           if (mpirank == mpisize - 1){ 
             int remainder = ntriangle%mpisize;
         if(remainder > 0) chunk += remainder;
           }
           end_triangle = chunk;
           char* data = new char[sizeof(double)*chunk*3*6];
           init_part_size = chunk;
           input_file.read(data, sizeof(double)*chunk*3*6);
           triangles = (double*) data;
           size_t size = chunk*nnpt;
           xt.resize(size);
           yt.resize(size);
           zt.resize(size);
           for(size_t i = 0; i < chunk; ++i) {
             for(size_t j = 0; j < nnpt; ++j){
               xt[i*nnpt+j] = triangles[i*nnpt*3+j*3+0];
               yt[i*nnpt+j] = triangles[i*nnpt*3+j*3+1];
               zt[i*nnpt+j] = triangles[i*nnpt*3+j*3+2];
             }
           }
         } else {
          std::cout<<"cannot read higher order mesh now" << std::endl;
        }
        input_file.close();       
      }
#else
      /**
      * A private unpartitioned reader for binary geometry data in I-DEAs format. 
      * Only rank 0 reads and broadcasts to the rest to supoort legacy behaviour
      */          
      void readBinaryIDEAsGeom(const char* file_name) { 
        if(mpirank == 0) {       
          std::ifstream input_file(file_name, std::ios::in|std::ios::binary|std::ios::ate);
          if(input_file.fail()) { 
            std::cout << "failed to open file " << file_name  <<  " by rank " << mpirank << " make sure file exists and access permissions are ok" <<std::endl;
            exit(EXIT_FAILURE);
         }
          std::stringstream ss;
          std::string line;
          int npat;
          double rcent[3];
          if(nnpt == 6) {
            char* firstRead = new char[sizeof(int)];
            input_file.seekg (0, std::ios::beg);
            input_file.read(firstRead, sizeof(int));
            nnod = *((int*)firstRead);
            input_file.read(firstRead, sizeof(int));
            npat = *((int*)firstRead);
            ntriangle=npat;
            init_part_size = ntriangle;
            sunod_ = new double[dim*nnod];
            nsupan_ = new int[npat*nnpt];
            sunod = new double*[dim];
            nsupan = new int*[npat];
            for (int i = 0; i < dim*nnod; ++i)
              sunod_[i] = 0.0;
            for (int i = 0; i < npat*nnpt; ++i)
              nsupan_[i] = 0;
            for (int i = 0; i < dim; ++i) 
              sunod[i] = (sunod_+i*nnod);        
            for (int i = 0; i < npat; ++i) 
              nsupan[i] = (nsupan_+i*nnpt);
            double sum1 = 0.0 , sum2 = 0.0, sum3 = 0.0;
            char* secondRead = new char[nnod*dim*sizeof(double)];
            input_file.read(secondRead, nnod*dim*sizeof(double));
            double* doubleValues = (double*)secondRead;//reinterpret as doubles

            for (int i = 0; i < nnod; ++i) {            
              sunod[0][i] = doubleValues[i*dim+0];
              sunod[1][i] = doubleValues[i*dim+1];
              sunod[2][i] = doubleValues[i*dim+2];
              sum1 += sunod[0][i];
              sum2 += sunod[1][i];
              sum3 += sunod[2][i];
            }
            rcent[0]=sum1/(double)nnod;
            rcent[1]=sum2/(double)nnod;
            rcent[2]=sum3/(double)nnod;
            for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < nnod; ++j) {              
                sunod[i][j] -= rcent[i];
              }
            }          
            char* thirdRead = new char[npat*nnpt*sizeof(int)]; 
            input_file.read(thirdRead,npat*nnpt*sizeof(int));
            int* intValues = (int*)thirdRead;
            for (int i = 0; i < npat; ++i) {
              nsupan[i][0] = intValues[i*nnpt+0];
              nsupan[i][5] = intValues[i*nnpt+1];
              nsupan[i][2] = intValues[i*nnpt+2];
              nsupan[i][4] = intValues[i*nnpt+3];
              nsupan[i][1] = intValues[i*nnpt+4];
              nsupan[i][3] = intValues[i*nnpt+5];
              for (int j = 0; j < nnpt; ++j) --nsupan[i][j];            
            }
            delete firstRead;
            delete secondRead;
            delete thirdRead;
        } else {
          std::cout<<"cannot read higher order mesh now" << std::endl;
        }       
        input_file.close();
      }
      sendRecvGeom();
    }
#endif
    /**
    * A private unpartitioned reader for textual geometry data in I-DEAs format. 
    * Only rank 0 reads and broadcasts to the rest to supoort legacy behaviour
    */              
    void readTextualIDEAsGeom(const char* file_name) { 
      if(mpirank == 0) {       
        std::ifstream input_file(file_name);
        if(input_file.fail()) { 
          std::cout << "failed to open file " << file_name  <<  " by rank " << mpirank << " make sure file exists and access permissions are ok" <<std::endl;
        }
        std::stringstream ss;
        std::string line;
        int npat;
        double rcent[3];
        if(nnpt == 6) {
          input_file >> nnod >> npat;
          std::getline(input_file,line);
          ntriangle=npat;
          init_part_size = ntriangle;
          sunod_ = new double[dim*nnod];
          nsupan_ = new int[npat*nnpt];
          sunod = new double*[dim];
          nsupan = new int*[npat];
          for (int i = 0; i < dim*nnod; ++i)
            sunod_[i] = 0.0;
          for (int i = 0; i < npat*nnpt; ++i)
            nsupan_[i] = 0;
          for (int i = 0; i < dim; ++i) 
            sunod[i] = (sunod_+i*nnod);        
          for (int i = 0; i < npat; ++i) 
            nsupan[i] = (nsupan_+i*nnpt);
          double sum1 = 0.0 , sum2 = 0.0, sum3 = 0.0;
          for (int i = 0; i < nnod; ++i) {            
            input_file >> sunod[0][i] >> sunod[1][i] >> sunod[2][i];                  
            sum1 += sunod[0][i];
            sum2 += sunod[1][i];
            sum3 += sunod[2][i];
          }
          
          rcent[0]=sum1/(double)nnod;
          rcent[1]=sum2/(double)nnod;
          rcent[2]=sum3/(double)nnod;
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < nnod; ++j) {              
              sunod[i][j] -= rcent[i];
            }
          }          
          for (int i = 0; i < npat; ++i) {
            input_file >> nsupan[i][0] >> nsupan[i][5]  >> nsupan[i][2] >> nsupan[i][4] >> nsupan[i][1] >> nsupan[i][3];
            for (int j = 0; j < nnpt; ++j) --nsupan[i][j];            
          }
      } else {
        std::cout<<"cannot read higher order mesh now" << std::endl;
      }       
      input_file.close();
    }
    sendRecvGeom();

  }
  /**
   * A private arguments reader for the command line args 
   * and the infrequently changable text args   
   */          
  void readInputArgs(const char* file_name) {
    std::ifstream input_file(file_name);
    if(input_file.fail()) { 
      std::cout << "failed to open file " << file_name  <<  " by rank " << mpirank << " make sure file exists and access permissions are ok" <<std::endl;
      exit(EXIT_FAILURE);
    }
    std::string line;
    double cd = 340.29;                          //sound speed of free space
    double the, phi, ethetad, ephid;
    double epsilon_r, mu_r; 
    double alpha,loss_sigma; 
    int order_tri;
    int order_basis;
    int nusolver;

    
    frequency = common_data::args->frequency;
    input_file >> the >> phi >> ethetad >> ephid;
    std::getline(input_file,line);
    
    the *= M_PI / 180.0;
    phi *= M_PI / 180.0;
    vkk[0] = std::sin(the) * std::cos(phi);
    vkk[1] = std::sin(the) * std::sin(phi);
    vkk[2] = std::cos(the);
    vee[0] = std::cos(the) * std::cos(phi) * ethetad - std::sin(phi) * ephid;
    vee[1] = std::cos(the) * std::sin(phi) * ethetad + std::cos(phi) * ephid;
    vee[2] =-std::sin(the) * ethetad;
    utils::cross_product(vkk,vee,vhh); 
    input_file >> the_sta_degree >> the_end_degree;
    std::getline(input_file,line);
    
    input_file >> phi_sta_degree >> phi_end_degree;
    std::getline(input_file,line);
    
    input_file >> nthe_rcs >> nphi_rcs;
    std::getline(input_file,line);
    
    input_file >> epsilon_r >> mu_r >> loss_sigma;
    std::getline(input_file,line);
    
    cd=cd / sqrt(epsilon_r * mu_r);            //light speed of background medium
    double wavenumkk = 2.0*M_PI*frequency/cd; //wave number k
    cjvk=wavenumkk;                            //j*k
    fmmAttributes.kreal = std::real(cjvk);
    fmmAttributes.kimag = std::imag(cjvk);  
    input_file >> alpha;
    std::getline(input_file,line);
    input_file >> order_tri;
    std::getline(input_file,line);

    if(order_tri == 1) 
      nnpt=3;
    else if(order_tri == 2) 
      nnpt=6;
    else
      std::cout<<"can not set higher order mesh now"<<std::endl;

    input_file >> order_basis;
    std::getline(input_file,line);

    if(order_basis == 0) 
      nipp = 1;
    else if(order_basis == 1) 
      nipp = 3;
    else if(order_basis == 2) 
      nipp = 6;
    else if(order_basis == 3) 
      nipp = 12;
    else
      std::cout<<"can not set higher order basis function now"<<std::endl;

    input_file >> nhdgqp;
    std::getline(input_file,line);

    input_file >> nlqp;
    std::getline(input_file,line);

    input_file >> nclgp;
    std::getline(input_file,line);    

    input_file >> nearpd;
    std::getline(input_file,line);

    input_file >> ngprcs;
    std::getline(input_file,line);

    input_file >> geomFormat;
    std::getline(input_file,line);

    input_file >> rnfobp;
    std::getline(input_file,line);

    input_file >> the_sta_degree_nf >> the_end_degree_nf;
    std::getline(input_file,line);

    input_file >> phi_sta_degree_nf >> phi_end_degree_nf;
    std::getline(input_file,line);

    input_file >> nthe_nf >> nphi_nf;
    std::getline(input_file,line);
    nitermax = common_data::args->maxiter;    
    precis = common_data::args->precision;
    fmmVerbose = common_data::args->fmmverbose;
    fmmAttributes.nthreads = common_data::args->threads;
    fmmAttributes.listbased = common_data::args->listbased;
    fmmAttributes.tspawn = common_data::args->nspawn;
    fmmAttributes.ncrit = common_data::args->ncrit;
    fmmAttributes.partitioning = common_data::args->partitioning;
    checkFMMDirect = common_data::args->direct;
    gmresRestart = common_data::args->gmresrestart;
    writeTimingOutputs = common_data::args->writeoutput;
    geomatryFileName = common_data::args->geomfile;

    input_file.close();
  }
};

struct integral_data
{ 
public:
  d_vector alphas;
  d_vector alphas1;
  d_vector alphas2;
  d_vector alphas3;
  d_vector betas;
  d_vector betas1; 
  d_vector betas2; 
  d_vector betas3; 
  d_vector gammas;
  d_vector gammas1;
  d_vector gammas2;
  d_vector gammas3;
  d_vector ws;
  d_vector ws1; 
  d_vector ws2; 
  d_vector ws3; 
  d_vector alphao;
  d_vector alphao1;
  d_vector alphao2;
  d_vector alphao3;
  d_vector betao;
  d_vector betao1; 
  d_vector betao2;
  d_vector betao3;
  d_vector gammao;
  d_vector gammao1;
  d_vector gammao2;
  d_vector gammao3;
  d_vector wo;
  d_vector wo1;
  d_vector wo2;
  d_vector wo3;
  d_vector wwwl;
  d_vector xxxl;
  d_vector wwo;
  d_vector xxo;
  d_vector wws;
  d_vector xxs;    
  int nipp, nhdgqp, ngprcs, gaussintpara1, gaussintpara2, nlqp, nclgp;
  int ngo, ngs, ngo1, ngs1, ngo2, ngs2, ngo3, ngs3, nlg, nlgo, nlgs;
  
  static integral_data& initIntegrationData(common_data const& comm_data) { 
    static integral_data instance(comm_data.nipp, comm_data.nhdgqp, comm_data.ngprcs,12,6,
      comm_data.nlqp, comm_data.nclgp);
    return instance;
  }

private:
  integral_data(int _nipp, int _nhdgqp, int _ngprcs, int _gaussintpara1, int _gaussintpara2, int _nlqp, int _nclgp) :
  nipp(_nipp), nhdgqp(_nhdgqp), ngprcs(_ngprcs),  gaussintpara1(_gaussintpara1), gaussintpara2(_gaussintpara2),
  nlqp(_nlqp), nclgp(_nclgp),
  ngo(0),  ngs(0), ngo1(0), ngs1(0), ngo2(0), ngs2(0), ngo3(0), 
  ngs3(0), nlg(0), nlgo(0), nlgs(0) {
    set_gauss_int_para(alphas,  betas,  gammas,  ws, alphao, betao, gammao, wo, ngo, ngs, nipp, nipp);
    set_gauss_int_para(alphas1, betas1, gammas1, ws1,alphao1,betao1,gammao1,wo1,ngo1,ngs1,nipp, nhdgqp);
    set_gauss_int_line(nlqp,nlgo,wwo,xxo);
    set_gauss_int_line(nlqp,nlgs,wws,xxs);
    set_gauss_int_para(alphas2,  betas2,  gammas2,  ws2, alphao2,  betao2,  gammao2,  wo2, ngo2, ngs2, nipp, ngprcs);
    set_gauss_int_line(nclgp,nlg,wwwl,xxxl);
    set_gauss_int_para(alphas3,  betas3,  gammas3,  ws3, alphao3, betao3, gammao3, wo3, ngo3, ngs3, gaussintpara1, gaussintpara2);
  }
  integral_data(integral_data const&);              
  void operator=(integral_data const&); 

  void ginttrg(int _n, d_vector& alpha, d_vector& beta, d_vector& gamma, d_vector& ww) {
    double v1, v2, v3, v4, v5, v7, v8, v9;
      if(_n == 1) { // p =1
       ww[0] = 1.0;
       v1 = 1.0/3.0;
       alpha[0] = v1; beta[0] = v1; gamma[0] = v1;
     }

      if(_n == 3) { // p =2
       std::fill(ww.begin(), ww.end() ,1.0/3.0);
       v1 = 2.0/3.0;
       v2 = 1.0/6.0;
       alpha[0] = v1; alpha[1] = v2;  alpha[2] = v2;  
       beta[0] = v2;  beta[1] = v1;   beta[2] = v2;    
       gamma[0] = v2; gamma[1] = v2;  gamma[2] = v1;   
     }

      if(_n == 4) { // p = 3
       v1 = 1.0/3.0;
       ww[0] = -0.5625;
       alpha[0] = v1; 
       beta[0] = v1;  
       gamma[0] = v1;
       std::fill(ww.begin()+1, ww.end() ,0.52083333333333333);
       alpha[1] = 0.6; beta[1] = 0.2;   gamma[1] = 0.2;  
       alpha[2] = 0.2; beta[2] = 0.6;   gamma[2] = 0.2;   
       alpha[3] = 0.2; beta[3] = 0.2;   gamma[3] = 0.6;
     }
      if(_n == 6) { // p =4
       std::fill(ww.begin()  , ww.begin()+3, 0.223381589678011);
       std::fill(ww.begin()+3, ww.end()    , 0.109951743655322);
       v1 = 0.108103018168070;
       v3 = 0.816847572980459;
       v2 = 0.445948490915965;
       v4 = 0.091576213509771;
       alpha[0]  =  v1;      beta[0]  =  v2;   gamma[0]  = v2;
       alpha[1]  =  v2;      beta[1]  =  v1;   gamma[1]  = v2;     
       alpha[2]  =  v2;      beta[2]  =  v2;   gamma[2]  = v1;  
       alpha[3]  =  v3;      beta[3]  =  v4;   gamma[3]  = v4;
       alpha[4]  =  v4;      beta[4]  =  v3;   gamma[4]  = v4;        
       alpha[5]  =  v4;      beta[5]  =  v4;   gamma[5]  = v3;    
     }
      if(_n == 7) { // p =5  //see in graglia paper
       ww[0]    = 0.225;
       alpha[0] = 1.0/3.0;
       beta[0]  = 1.0/3.0;
       gamma[0] = 1.0/3.0;
       std::fill(ww.begin()+1  , ww.begin()+4, (155.0+sqrt(15.0))/1200.0);
       std::fill(ww.begin()+4  , ww.end(),     (155.0-sqrt(15.0))/1200.0);
       v1 = (9.0-2.0*sqrt(15.0))/21.0;
       v2 = (6.0-   sqrt(15.0))/21.0;
       v3 = (6.0+   sqrt(15.0))/21.0;
       v4 = (9.0+2.0*sqrt(15.0))/21.0;
       alpha[1] = v1;  alpha[2] = v3;  alpha[3] = v3; 
       beta[1] = v3;   beta[2] = v1;   beta[3] = v3; 
       gamma[1] = v3;  gamma[2] = v3;  gamma[3] = v1;
       alpha[4] = v4;  alpha[5] = v2;  alpha[6] = v2;
       beta[4] = v2;   beta[5] = v4;   beta[6] = v2;
       gamma[4] = v2;  gamma[5] = v2;  gamma[6] = v4;
     }
      if(_n == 9) { // p =3  
       std::fill(ww.begin()  , ww.end(), 1.0/9.0);
       v1 =  1.0/9.0;
       v2 =  2.0/9.0;
       v3 =  4.0/9.0;
       v4 =  5.0/9.0;
       v5 =  7.0/9.0;
       alpha[0] = v1;  beta[0] = v5; gamma[0] = 1.0 - beta[0] - alpha[0];
       alpha[1] = v5;  beta[1] = v1; gamma[1] = 1.0 - beta[1] - alpha[1]; 
       alpha[2] = v1;  beta[2] = v1; gamma[2] = 1.0 - beta[2] - alpha[2];
       alpha[3] = v3;  beta[3] = v1; gamma[3] = 1.0 - beta[3] - alpha[3]; 
       alpha[4] = v1;  beta[4] = v3; gamma[4] = 1.0 - beta[4] - alpha[4];
       alpha[5] = v3;  beta[5] = v3; gamma[5] = 1.0 - beta[5] - alpha[5];
       alpha[6] = v2;  beta[6] = v4; gamma[6] = 1.0 - beta[6] - alpha[6];
       alpha[7] = v2;  beta[7] = v2; gamma[7] = 1.0 - beta[7] - alpha[7];
       alpha[8] = v4;  beta[8] = v2; gamma[8] = 1.0 - beta[8] - alpha[8];
     }

      if(_n == 12) { // p = 6
       std::fill(ww.begin(),   ww.begin()+3,  0.116786275726379);
       std::fill(ww.begin()+3, ww.begin()+6,  0.050844906370207);
       std::fill(ww.begin()+6, ww.end(),      0.082851075618374);

       v1  = 0.501426509658179;  v2 = 0.249286745170910;
       v4  = 0.873821971016996;  v5 = 0.063089014491502;
       v7  = 0.053145049844817;  v8 = 0.310352451033784;
       v9  = 1.0 - v7 - v8;

       alpha[0]  =  v1;  beta[0]  =  v2;    gamma[0]  =  v2;
       alpha[1]  =  v2;  beta[1]  =  v1;    gamma[1]  =  v2;       
       alpha[2]  =  v2;  beta[2]  =  v2;    gamma[2]  =  v1;     
       alpha[3]  =  v4;  beta[3]  =  v5;    gamma[3]  =  v5;
       alpha[4]  =  v5;  beta[4]  =  v4;    gamma[4]  =  v5;         
       alpha[5]  =  v5;  beta[5]  =  v5;    gamma[5]  =  v4;        
       alpha[6]  =  v7;  beta[6]  =  v8;    gamma[6]  =  v9;
       alpha[7]  =  v7;  beta[7]  =  v9;    gamma[7]  =  v8;    
       alpha[8]  =  v8;  beta[8]  =  v7;    gamma[8]  =  v9;      
       alpha[9]  =  v8;  beta[9]  =  v9;    gamma[9]  =  v7;      
       alpha[10] =  v9;  beta[10]  =  v8;   gamma[10]  =  v7;     
       alpha[11] =  v9;  beta[11]  =  v7;   gamma[11]  =  v8;    
     }

      if(_n == 13) { // p =7
       ww[0]     =   -0.149570044467682;
       std::fill(ww.begin()+1, ww.begin()+4, 0.175615257433208);
       std::fill(ww.begin()+4, ww.begin()+7, 0.053347235608838);
       std::fill(ww.begin()+7, ww.end(),     0.077113760890257);

       alpha[0]  = 0.333333333333333;  beta[0] = 0.333333333333333;
       alpha[1]  = 0.479308067841920;  beta[1] = 0.260345966079040;
       alpha[4]  = 0.869739794195568;  beta[4] = 0.065130102902216;
       alpha[7]  = 0.048690315425316;  beta[7] = 0.312865496004874;

       gamma[0]  = 1.0 - alpha[0] - beta[0];
       gamma[1]  = 1.0 - alpha[1] - beta[1];
       gamma[4]  = 1.0 - alpha[4] - beta[4];
       gamma[7]  = 1.0 - alpha[7] - beta[7];

       alpha[2]  =  beta[1];         beta[2]  = alpha[1];           
       alpha[3]  =  beta[1];         beta[3]  =  beta[1];           
       alpha[5]  =  beta[4];         beta[5]  = alpha[4];           
       alpha[6]  =  beta[4];         beta[6]  =  beta[4];           

       alpha[8]   = alpha[7];         beta[8]  = gamma[7];           
       alpha[9]   =  beta[7];         beta[9]  = alpha[7];           
       alpha[10]  =  beta[7];        beta[10]  = gamma[7];           
       alpha[11]  = gamma[7];        beta[11]  =  beta[7];           
       alpha[12]  = gamma[7];        beta[12]  = alpha[7];

       gamma[2]     = 1.0 - alpha[2] - beta[2];
       gamma[3]     = 1.0 - alpha[3] - beta[3];
       gamma[5]     = 1.0 - alpha[5] - beta[5];
       gamma[6]     = 1.0 - alpha[6] - beta[6];

       gamma[8]     = 1.0 - alpha[8] - beta[8];
       gamma[9]     = 1.0 - alpha[9] - beta[9];
       gamma[10]    = 1.0 - alpha[10] - beta[10];
       gamma[11]    = 1.0 - alpha[11] - beta[11];
       gamma[12]    = 1.0 - alpha[12] - beta[12];
     }

      if(_n == 16) { // p = 8
        std::fill(ww.begin(),   ww.begin()+3,  0.095091634267285);
        std::fill(ww.begin()+3, ww.begin()+6,  0.103217370534718);
        std::fill(ww.begin()+6, ww.begin()+12, 0.027230314174435);
        std::fill(ww.begin()+12,ww.begin()+15, 0.032458497623198);
        ww[15]= 0.144315607677787;
        v1  = 0.081414823414554;  v2 = 0.459292588292723;
        v4  = 0.658861384496480;  v5 = 0.170569307751760;
        v7  = 0.008394777409958;  v8 = 0.263112829634638;
        v9  = 1.0 - v7 - v8;
        alpha[ 0]  =  v1;  beta[ 0]  =  v2;    gamma[ 0]  =  v2;
        alpha[ 1]  =  v2;  beta[ 1]  =  v1;    gamma[ 1]  =  v2;       
        alpha[ 2]  =  v2;  beta[ 2]  =  v2;    gamma[ 2]  =  v1;     
        alpha[ 3]  =  v4;  beta[ 3]  =  v5;    gamma[ 3]  =  v5;
        alpha[ 4]  =  v5;  beta[ 4]  =  v4;    gamma[ 4]  =  v5;         
        alpha[ 5]  =  v5;  beta[ 5]  =  v5;    gamma[ 5]  =  v4;        
        alpha[ 6]  =  v7;  beta[ 6]  =  v8;    gamma[ 6]  =  v9;
        alpha[ 7]  =  v7;  beta[ 7]  =  v9;    gamma[ 7]  =  v8;    
        alpha[ 8]  =  v8;  beta[ 8]  =  v7;    gamma[ 8]  =  v9;      
        alpha[ 9]  =  v8;  beta[ 9]  =  v9;    gamma[ 9]  =  v7;      
        alpha[10]  =  v9;  beta[10]  =  v8;    gamma[10]  =  v7;     
        alpha[11]  =  v9;  beta[11]  =  v7;    gamma[11]  =  v8;    
        v1  = 0.898905543365938;  v2 = 0.050547228317031;
        alpha[12]  =  v1;  beta[12]  =  v2;    gamma[12]  =  v2;
        alpha[13]  =  v2;  beta[13]  =  v1;    gamma[13]  =  v2;       
        alpha[14]  =  v2;  beta[14]  =  v2;    gamma[14]  =  v1;   
        v1 = 1.0/3.0;
        alpha[15]  =  v1;  beta[15]  =  v1;    gamma[15]  =  v1;   
      }
      
      if(_n == 25) { // p = 10
        std::fill(ww.begin(),   ww.begin()+3,  0.036725957756467);
        std::fill(ww.begin()+3, ww.begin()+6,  0.045321059435528);
        std::fill(ww.begin()+6, ww.begin()+12, 0.072757916845420);
        std::fill(ww.begin()+12,ww.begin()+18, 0.028327242531057);
        std::fill(ww.begin()+18,ww.begin()+24, 0.009421666963733);
        ww[24] = 0.090817990382754;

        v1  = 0.028844733232685;  v2 = 0.485577633383657;
        v4  = 0.781036849029926;  v5 = 0.109481575485037;
        v7  = 0.141707219414880;  v8 = 0.307939838764121;
        v9  = 1.0 - v7 - v8;

        alpha[ 0]  =  v1;  beta[ 0]  =  v2;    gamma[ 0]  =  v2;
        alpha[ 1]  =  v2;  beta[ 1]  =  v1;    gamma[ 1]  =  v2;       
        alpha[ 2]  =  v2;  beta[ 2]  =  v2;    gamma[ 2]  =  v1;     
        alpha[ 3]  =  v4;  beta[ 3]  =  v5;    gamma[ 3]  =  v5;
        alpha[ 4]  =  v5;  beta[ 4]  =  v4;    gamma[ 4]  =  v5;         
        alpha[ 5]  =  v5;  beta[ 5]  =  v5;    gamma[ 5]  =  v4;        
        alpha[ 6]  =  v7;  beta[ 6]  =  v8;    gamma[ 6]  =  v9;
        alpha[ 7]  =  v7;  beta[ 7]  =  v9;    gamma[ 7]  =  v8;    
        alpha[ 8]  =  v8;  beta[ 8]  =  v7;    gamma[ 8]  =  v9;      
        alpha[9]  =  v8;  beta[9]  =  v9;    gamma[9]  =  v7;      
        alpha[10] =  v9;  beta[10]  =  v8;    gamma[10]  =  v7;     
        alpha[11] =  v9;  beta[11]  =  v7;    gamma[11]  =  v8;    
        v7  = 0.025003534762686;  v8 = 0.246672560639903;
        v9  = 1.0 - v7 - v8;
        alpha[12]  =  v7;  beta[12]  =  v8;    gamma[12]  =  v9;
        alpha[13]  =  v7;  beta[13]  =  v9;    gamma[13]  =  v8;    
        alpha[14]  =  v8;  beta[14]  =  v7;    gamma[14]  =  v9;      
        alpha[15]  =  v8;  beta[15]  =  v9;    gamma[15]  =  v7;      
        alpha[16]  =  v9;  beta[16]  =  v8;    gamma[16]  =  v7;     
        alpha[17]  =  v9;  beta[17]  =  v7;    gamma[17]  =  v8;      
        v7  = 0.009540815400299;  v8 = 0.066803251012200;
        v9  = 1.0 - v7 - v8;
        alpha[18]  =  v7;  beta[18]  =  v8;    gamma[18]  =  v9;
        alpha[19]  =  v7;  beta[19]  =  v9;    gamma[19]  =  v8;    
        alpha[20]  =  v8;  beta[20]  =  v7;    gamma[20]  =  v9;      
        alpha[21]  =  v8;  beta[21]  =  v9;    gamma[21]  =  v7;      
        alpha[22]  =  v9;  beta[22]  =  v8;    gamma[22]  =  v7;     
        alpha[23]  =  v9;  beta[23]  =  v7;    gamma[23]  =  v8;    
        v1 = 1.0/3.0;
        alpha[24]  =  v1;  beta[24]  =  v1;    gamma[24]  =  v1;   
      }
      if(_n == 42) { // p =14 
        std::fill(ww.begin(),   ww.begin()+3,  0.021883581369429);
        std::fill(ww.begin()+3, ww.begin()+6,  0.032788353544125);
        std::fill(ww.begin()+6, ww.begin()+9,  0.051774104507292);
        std::fill(ww.begin()+9, ww.begin()+12, 0.042162588736993);
        std::fill(ww.begin()+12,ww.begin()+15, 0.014433699669777);
        std::fill(ww.begin()+15,ww.begin()+18, 0.004923403602400);
        std::fill(ww.begin()+18,ww.begin()+24, 0.024665753212564);
        std::fill(ww.begin()+24,ww.begin()+30, 0.038571510787061);
        std::fill(ww.begin()+30,ww.begin()+36, 0.014436308113534);
        std::fill(ww.begin()+36,  ww.end(),    0.005010228838501);

        alpha[ 0] = 0.022072179275643;  beta[ 0] = (1.0 - alpha[0])/2.0; 
        alpha[ 3] = 0.164710561319092;  beta[ 3] = (1.0 - alpha[3])/2.0;
        alpha[ 6] = 0.453044943382323;  beta[ 6] = (1.0 - alpha[6])/2.0;
        alpha[ 9] = 0.645588935174913;  beta[ 9] = (1.0 - alpha[9])/2.0;
        alpha[12] = 0.876400233818255;  beta[12] = (1.0 - alpha[12])/2.0;
        alpha[15] = 0.961218077502598;  beta[15] = (1.0 - alpha[15])/2.0;


        alpha[18] = 0.057124757403648;  beta[18] = 0.172266687821356;
        alpha[24] = 0.092916249356972;  beta[24] = 0.336861459796345;
        alpha[30] = 0.014646950055654;  beta[30] = 0.298372882136258;
        alpha[36] = 0.001268330932872;  beta[36] = 0.118974497696957;

        for (int i = 0; i < _n; ++i) 
          gamma[i]  = 1.0 - alpha[i] - beta[i];

        alpha[ 1]  =  beta[ 0];         beta[ 1]  = alpha[ 0];           
        alpha[ 2]  =  beta[ 0];         beta[ 2]  = gamma[ 0];

        alpha[ 4]  =  beta[ 3];         beta[ 4]  = alpha[ 3];          
        alpha[ 5]  =  beta[ 3];         beta[ 5]  = gamma[ 3]; 

        alpha[ 7]  =  beta[ 6];         beta[ 7]  = alpha[ 6];          
        alpha[ 8]  =  beta[ 6];         beta[ 8]  = gamma[ 6];

        alpha[10]  =  beta[ 9];         beta[10]  = alpha[ 9];          
        alpha[11]  =  beta[ 9];         beta[11]  = gamma[ 9];  

        alpha[13]  =  beta[12];         beta[13]  = alpha[12]; 
        alpha[14]  =  beta[12];         beta[14]  = gamma[12]; 

        alpha[16]  =  beta[15];         beta[16]  = alpha[15];          
        alpha[17]  =  beta[15];         beta[17]  = gamma[15]; 

        alpha[19]  = alpha[18];         beta[19]  = gamma[18];          
        alpha[20]  =  beta[18];         beta[20]  = alpha[18];          
        alpha[21]  =  beta[18];         beta[21]  = gamma[18];          
        alpha[22]  = gamma[18];         beta[22]  =  beta[18];          
        alpha[23]  = gamma[18];         beta[23]  = alpha[18]; 

        alpha[25]  = alpha[24];         beta[25]  = gamma[24];          
        alpha[26]  =  beta[24];         beta[26]  = alpha[24];          
        alpha[27]  =  beta[24];         beta[27]  = gamma[24];          
        alpha[28]  = gamma[24];         beta[28]  =  beta[24];          
        alpha[29]  = gamma[24];         beta[29]  = alpha[24];  

        alpha[31]  = alpha[30];         beta[31]  = gamma[30];          
        alpha[32]  =  beta[30];         beta[32]  = alpha[30];          
        alpha[33]  =  beta[30];         beta[33]  = gamma[30];          
        alpha[34]  = gamma[30];         beta[34]  =  beta[30];          
        alpha[35]  = gamma[30];         beta[35]  = alpha[30]; 

        alpha[37]  = alpha[36];         beta[37]  = gamma[36];          
        alpha[38]  =  beta[36];         beta[38]  = alpha[36];          
        alpha[39]  =  beta[36];         beta[39]  = gamma[36];          
        alpha[40]  = gamma[36];         beta[40]  =  beta[36];          
        alpha[41]  = gamma[36];         beta[41]  = alpha[36];   
        
        for (int i = 0; i < _n; ++i) 
          gamma[i]  = 1.0 - alpha[i] - beta[i];
      }

    }

    void set_gauss_int_para(d_vector& _alphas,  d_vector& _betas,  d_vector& _gammas,  d_vector& _ws, 
     d_vector& _alphao, d_vector& _betao, d_vector& _gammao, d_vector& _wo,
     int& _ngo, int& _ngs, int _ng1, int _ng2) {
      if(_wo.size() > 0 && _ngo == _ng1 && _ngs == _ng2 ) return;
      _ngo = _ng1;
      _ngs = _ng2;
      assert(_ngo > 0);
      assert(_ngs > 0);
      _alphas.resize(_ngs,0);
      _betas.resize(_ngs,0);
      _gammas.resize(_ngs,0);
      _ws.resize(_ngs,0);
      _alphao.resize(_ngo,0);
      _betao.resize(_ngo,0);
      _gammao.resize(_ngo,0);
      _wo.resize(_ngo,0);
      ginttrg(_ngo, _alphao, _betao, _gammao, _wo);
      ginttrg(_ngs, _alphas, _betas, _gammas, _ws);
    }

    void gau_leg(double x1, double x2, d_vector& x, d_vector& w, int n) {
      double const eps = 3.0e-14;
      double p1,p2,p3,pp,z,z1;
      int m = (n+1)/2;
      double xm=0.5*(x2+x1);
      double xl=0.5*(x2-x1);
      for (int i = 1; i <= m; ++i) {
        z=cos(3.141592654*(i-.25)/(n+.5));
        do {
          p1=1.0;
          p2=0.0;
          for (int j = 1; j <= n; ++j) {
            p3=p2;
            p2=p1;
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
          }
          pp=n*(z*p1-p2)/(z*z-1.0);
          z1=z;
          z=z1-p1/pp; 
        } while(std::abs(z-z1) > eps);
        x[i-1]=xm-xl*z;
        x[n-i]=xm+xl*z;
        w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n-i]=w[i-1];
      }
    }

    void set_gauss_int_line(int _nlgin, int& _nlg, d_vector& _www, d_vector& _xxx) {
      assert(_nlgin>0);
      if(_www.size() > 0 && _nlg > _nlgin) return;
      _nlg = _nlgin;
      _www.resize(_nlgin,0.0);
      _xxx.resize(_nlgin,0.0);
      gau_leg( 0.0,1.0,_xxx,_www,_nlgin);
    }
  };
  Args* common_data::args = NULL;
} // end bemfmm workspace

#endif
