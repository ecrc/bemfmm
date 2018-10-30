.SUFFIXES: .o .c .C .h .cxx .cpp
.PHONY: clean intel shaheen_gnu shaheen_intel gnu

BLAS_PATH=/home/abduljm/libs/petsc-3.9.3/arch-linux2-c-debug
BLAS_LIB=${BLAS_PATH}/lib

LAPACK_PATH=/home/abduljm/libs/petsc-3.9.3/arch-linux2-c-debug
LAPACK_LIB=${LAPACK_PATH}/lib

LIB_BLAS=-L${BLAS_LIB} -lfblas
LIB_LAPACK=-L${LAPACK_LIB} -lflapack 

LD_FLAGS+=${LIB_LAPACK} ${LIB_BLAS} -lgfortran -lm


# METIS
METIS_PATH = METIS
METIS_LIBS = -L${METIS_PATH}/lib/ -lparmetis -L${METIS_PATH}/lib/ -lmetis
METIS_INCS = -I${METIS_PATH}/inc

# TBB
TBB_PATH = TBB
TBB_LIBS = -L${TBB_PATH}/lib/ -ltbb
TBB_INCS = -I${TBB_PATH}/inc/

# FMM
P = 10
FMM_FLGS =  -DEXAFMM_EXPANSION=$P \
            -DEXAFMM_SPHERICAL \
            -DEXAFMM_HELMHOLTZ \
            -DEXAFMM_ACOUSTICS \
            -DEXAFMM_NEARF_TREATMENT \
            -DEXAFMM_COUNT_KERNEL \
            -DEXAFMM_USE_PARMETIS \
            -DEXAFMM_WITH_TBB \
            -DUSE_FMM #\
            -DUSE_PART

FMM_INCS = -Iinclude/fmm

# Core
LIBS = ${METIS_LIBS} ${TBB_LIBS} ${LD_FLAGS}
INCS = -Iinclude/ ${METIS_INCS} ${TBB_INCS} ${FMM_INCS}
FLGS = ${INCS} ${FMM_FLGS}  -mavx -O3
CXX = mpicxx

intel: CXX = mpiicpc
intel: FLGS += -qopenmp
intel: clean bemfmm_test_mpi

gnu: CXX = mpicxx
gnu: FLGS += -fopenmp
gnu: clean bemfmm_test_mpi

shaheen: CXX = CC
shaheen: LIBS += -dynamic
shaheen: FLGS += -fopenmp
shaheen: clean bemfmm_test_mpi

SRC = bemfmm_test_mpi.cxx
OBJ = $(SRC:.cxx=.o)

.cxx.o:
	$(CXX) $(FLGS) -c $<

bemfmm_test_mpi: $(OBJ)
	$(CXX) $(FLGS) $(OBJ) $(LIBS) -o $@ 

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.time
	rm -f *.mod
	rm -f *.o.optdbg
	rm -f *.bak
	rm -f *~	
	rm -f bemfmm_test_mpi
