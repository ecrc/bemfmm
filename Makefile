all: clean bemfmm_test_mpi

.SUFFIXES: .o .c .C .h .cxx .cpp
.PHONY: clean all

include make.inc

# METIS
METIS_LIBS = -L${METIS_LIB_PATH} -lparmetis -lmetis
METIS_INCS = -I${METIS_INC_PATH}

# TBB
TBB_LIBS = -L${TBB_LIB_PATH} -ltbb
TBB_INCS = -I${TBB_INC_PATH}

# LAPACK
LAPACK_LIBS = -L${LAPACK_LIB_PATH} -llapack -lrefblas -lgfortran -lm
BLAS_LIBS = -L${BLAS_LIB_PATH} -lrefblas

# Core
LIBS = ${METIS_LIBS} ${TBB_LIBS} ${LAPACK_LIBS} ${USERLIBS} ${BLAS_LIBS}
INCS = -Iinclude/fmm -Iinclude/ ${METIS_INCS} ${TBB_INCS} ${USERINCS}
FLGS = ${INCS} ${FMM_FLGS}

ifeq ($(MODE), DEV)
  FLGS += -g
else
  FLGS += -O3 -mavx
endif

ifeq ($(LINKING_TYPE), DYNAMIC)
	LIBS += -dynamic
endif

ifeq ($(ENV), GNU)
	FLGS += -fopenmp
else ifeq ($(ENV), INTEL)
	FLGS += -qopenmp
endif

FLGS += $(USERCXXFLAGS)

SRC = bemfmm_test_mpi.cxx
OBJ = $(SRC:.cxx=.o)

.cxx.o:
	$(CXX) $(FLGS) -c $<

bemfmm_test_mpi: $(OBJ)
	$(CXX) $(FLGS) $(OBJ) $(LIBS) -o $@
	$(RM) -f *.o

test_parallel:
	${EXEC} -n 4 ./bemfmm_test_mpi -f geom/sphere/geo_mesh_3156.inp -wd -t 10 -c 500 -q 683.13 -p h -r 30 -m 1000

test_serial:
	${EXEC} -n 1 ./bemfmm_test_mpi -f geom/sphere/geo_mesh_156.inp -wd -t 10 -c 500 -q 683.13 -p h -r 30 -m 1000 

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.time
	rm -f *.mod
	rm -f *.o.optdbg
	rm -f *.bak
	rm -f *~
	rm -f *.dat
	rm -f bemfmm_test_mpi
