/**
 *
 * @file args.h
 *
 * @copyright 2018 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 *
 * @author Mustafa Abduljabbar [mustafa.abduljabbar@kaust.edu.sa] and Mohammed Al Farhan [mohammed.farhan@kaust.edu.sa].
 *
 **/

#ifndef args_h
#define args_h
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <getopt.h>

extern int opterr;

namespace bemfmm {
static struct option long_options[] = {
  {"partitioning", required_argument, 0, 'p'},
  {"geomfile",     required_argument, 0, 'f'},
  {"configfile",   required_argument, 0, 'i'},
  {"threads",      required_argument, 0, 't'},
  {"frequency",    required_argument, 0, 'q'},
  {"precision",    required_argument, 0, 'e'},
  {"fmmverbose",   no_argument,       0, 'v'},
  {"writeoutput",  no_argument,       0, 'w'},
  {"ncrit",        required_argument, 0, 'c'}, 
  {"restart",      required_argument, 0, 'r'}, 
  {"checkdirect",  no_argument,       0, 'd'},
  {"listbased",    no_argument,       0, 'l'},
  {"nspawn",       required_argument, 0, 's'},
  {0, 0, 0, 0}
};
class Args {
public:  
  std::string geomfile;
  std::string configfile;
  std::string partitioning;
  int threads;
  int fmmverbose;
  int writeoutput;
  int ncrit;
  int maxiter;
  int gmresrestart; 
  int direct;
  int listbased;
  int nspawn;
  double frequency;
  double precision;

private:
  void usage(char * name) {
    fprintf(stderr,
            "Usage: %s [options]\n"
            "Long option (short option)     : Description (Default value)\n"
            " --partitioning (-p)[h/o/b/p]  : Hilbert, octsection, bisection, ParMETIS (%s)\n"
            " --geomfile (-f)               : The name of the file containing geometry (%s)\n"
            " --configfile (-i)             : The name of the file containing acoustics configration (%s)\n"
            " --threads (-t)                : Number of threads used in FMM traversal/tree building (%d)\n"
            " --frequency (-q)              : Single wave frequency  (%f)\n"            
            " --precision (-e)              : Iterative solver precision  (%f)\n"            
            " --fmmverbose (-v)             : Verbose FMM timing output (%d)\n"
            " --writeoutput (-w)            : Write results to files (%d)\n"
            " --ncrit (-c)                  : Number of bodies per leaf cell (FMM) (%d)\n"
            " --maxiter (-m)                : Max number of GMRES iterations (%d)\n"
            " --restart (-r)                : GMRES restart (%d)\n"
            " --direct (-d)                 : Verify FMM against direct  (%d)\n"
            " --listbased (-l)              : use list-based traversal  (%d)\n"            
            " --nspawn (-s)                 : Threshold for stopping task creation during recursion  (%d)\n"
            " --help (-h)                   : Show this help document\n",     
            name,
            partitioning.c_str(),
            geomfile.c_str(),
            configfile.c_str(),
            threads,
            frequency,
            precision,
            fmmverbose,
            writeoutput,
            ncrit, 
            maxiter,
            gmresrestart,
            direct, 
            listbased, 
            nspawn);
  }


public:
  Args(int argc = 0, char ** argv = NULL) :
    geomfile("geom/sphere/geo_mesh_156.inp"),
    configfile("application_params"),
    partitioning("b"),
    threads(4),
    fmmverbose(0),
    writeoutput(0),
    ncrit(500),
    maxiter(1000),
    gmresrestart(30),
    direct(0), 
    listbased(0), 
    nspawn(1000), 
    frequency(100.0),
    precision(1.0e-03){
    int cont = 1;
    opterr = 0;
    while (cont) {
      int option_index;
      int c = getopt_long(argc, argv, "s:m:r:f:e:p:i:t:q:vdwlc:h", long_options, &option_index);
      if (c == -1) break;
      switch (c) {
      case 'c':
        ncrit = atoi(optarg);
        break;
      case 'm':
        maxiter = atoi(optarg);
        break;        
      case 'r':
        gmresrestart = atoi(optarg);
        break;        
      case 'q':
        frequency = atof(optarg);
        break;
      case 'e':
        precision = atof(optarg);
        break;        
      case 't':
        threads = atoi(optarg);
        break; 
      case 's':
        nspawn = atoi(optarg);
        break;         
      case 'v':
        fmmverbose = 1;
        break;                
      case 'w':
        writeoutput = 1;
        break;        
      case 'd':
        direct = 1;
      break;  
      case 'l':
        listbased = 1;
      break;  
      case 'i':
        configfile = std::string(optarg);
        break;      
      case 'f':
        geomfile   = std::string(optarg);
      break;      
      case 'p':
        partitioning  = std::string(optarg);
      break;      
      case 'h':
        usage(argv[0]);
        exit(0);      
      default:
        cont = 0;
        break;
      }
    }
  }
};
}

#endif
