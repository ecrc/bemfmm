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
  {"fmmverbose",   no_argument,       0, 'v'},
  {"writeoutput",  no_argument,       0, 'w'},
  {"ncrit",        required_argument, 0, 'c'}, 
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
  int direct;
  int listbased;
  int nspawn;
  double frequency;

private:
  void usage(char * name) {
    fprintf(stderr,
            "Usage: %s [options]\n"
            "Long option (short option)     : Description (Default value)\n"
            " --partitioning (-f)           : The FMM partitioning used (%s)\n"
            " --geomfile (-f)               : The name of the file containing geometry (%s)\n"
            " --configfile (-i)             : The name of the file containing acoustics configration (%s)\n"
            " --threads (-t)                : Number of threads used in FMM traversal/tree building (%d)\n"
            " --frequency (-q)              : Single wave frequency  (%f)\n"            
            " --fmmverbose (-v)             : Verbose FMM output (%d)\n"
            " --writeoutput (-w)            : Write Acoustics output to files (%d)\n"
            " --ncrit (-c)                  : Number of bodies per leaf cell (FMM) (%d)\n"
            " --direct (-d)                 : run direct kernel using N^2  (%d)\n"
            " --listbased (-l)              : use list-based traversal  (%d)\n"            
            " --nspawn (-s)                 : Threshold for stopping task creation during recursion  (%d)\n"
            " --help (-h)                   : Show this help document\n",     
            name,
            partitioning.c_str(),
            geomfile.c_str(),
            configfile.c_str(),
            threads,
	    frequency,
            fmmverbose,
            writeoutput,
            ncrit, 
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
    direct(0), 
    listbased(0), 
    nspawn(1000), 
    frequency(100){
    int cont = 1;
    opterr = 0;
    while (cont) {
      int option_index;
      int c = getopt_long(argc, argv, "s:f:p:i:t:q:vdwlc:h", long_options, &option_index);
      if (c == -1) break;
      switch (c) {
      case 'c':
        ncrit = atoi(optarg);
        break;
      case 'q':
        frequency = atof(optarg);
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
