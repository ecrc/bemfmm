
all: clean bemfmm_test_mpi

.SUFFIXES: .o .c .C .h .cxx .cpp
.PHONY: clean intel shaheen_gnu shaheen_intel gnu

include make.inc

# METIS
METIS_LIBS = -L${METIS_LIB_PATH} -lparmetis -lmetis
METIS_INCS = -I${METIS_INC_PATH}

# TBB
TBB_LIBS = -L${TBB_LIB_PATH} -ltbb
TBB_INCS = -I${TBB_INC_PATH}

# LAPACK
LAPACK_LIBS = -L${LAPACK_LIB_PATH} -lflapack -L${LAPACK_LIB_PATH} -lfblas -lm -lgfortran

# Core
LIBS = ${METIS_LIBS} ${TBB_LIBS} ${LAPACK_LIBS}
INCS = -Iinclude/fmm -Iinclude/ ${METIS_INCS} ${TBB_INCS}
FLGS = ${INCS} ${FMM_FLGS}

ifeq ($(MODE), DEV)
  FLGS += -g  
else
  FLGS += -O3 -mavx
endif

ifeq ($(COMP), SHAHEEN)
	CXX = CC
	LIBS += -dynamic
	ifeq (${ENV}, SHAHEEN-GNU)
		FLGS += -fopenmp
	else ifeq (${ENV}, SHAHEEM-INTEL)
		FLGS += -qopenmp
	endif
else ifeq ($(COMP), GNU)
	CXX = mpicxx
	FLGS += -fopenmp
else ifeq ($(COMP), INTEL)
	CXX = mpicpc
	FLGS += -qopenmp
endif

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
	rm -f *.dat
	rm -f *.o.optdbg
	rm -f *.bak
	rm -f *~	
	rm -f bemfmm_test_mpi
