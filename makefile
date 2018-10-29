#===========================GNU========================================
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

ifeq ($(THREADING_MODEL),)
$(error THREADING_MODEL is not set)
endif

ifneq (,$(findstring tbb,$(THREADING_MODEL)))
ifeq ($(TBBROOT),)
$(error TBBROOT is not set)
endif
CXXFLAGS+=-DEXAFMM_WITH_TBB
#LD_FLAGS+=-L$TBBROOT/lib/intel64/gcc4.7/ -Wl,-R/opt/intel/composer_xe_2015.2.164/compiler/lib/intel64/ -lifcore -dynamic
LD_FLAGS+=-dynamic -ltbb

else ifneq (,$(findstring cilk,$(THREADING_MODEL)))
CXXFLAGS+=-DEXAFMM_WITH_CILK
else ifneq (,$(findstring mthread,$(THREADING_MODEL)))
ifeq ($(MTHREADROOT),)
$(error MTHREADROOT is not set)
endif
CXXFLAGS+=-DEXAFMM_WITH_MTHREAD -I$(MTHREADROOT)/include
LD_FLAGS+=-L$(MTHREADROOT)/lib -Wl,-R$(MTHREADROOT)/lib -lmyth
else
CXXFLAGS+=-DEXAFMM_WITH_OMP
endif 

ifeq ($(parmetis),1)
	PARMETIS_PATH = /home/abduljm/Repos/parmetis-4.0.3/
	LIB_PARMETIS = -L${PARMETIS_PATH}/build/Linux-x86_64/libparmetis/ -lparmetis  -L${PARMETIS_PATH}/build/Linux-x86_64/libmetis/ -lmetis
endif

#ifeq ($(EXAFMM_PATH),)
#$(error EXAFMM_PATH is not set)
#endif
ifeq ($(P),)
P = 10
$(info $$$$$$$$$$$$$$$$$$$$  P is not set, setting it to $P  $$$$$$$$$$$$$$$$$$$$$$$) 
endif

FMM_FLAGS=-DEXAFMM_EXPANSION=$P -DEXAFMM_SPHERICAL -DEXAFMM_HELMHOLTZ -DEXAFMM_ACOUSTICS -DEXAFMM_NEARF_TREATMENT -DEXAFMM_COUNT_KERNEL #-DEXAFMM_USE_PARMETIS #-DEXAFMM_USE_P2P# -DEXAFMM_NO_M2L #-DEXAFMM_SINGLE #-DEXAFMM_USE_AUTOVEC #-DEXAFMM_NO_M2L -DEXAFMM_SINGLE
CXXFLAGS+=${FMM_FLAGS}


ifeq ($(ACOUS_ENV),)
$(info $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ACOUS_ENV is not set. Set it to either gnu or intel. gnu will be used by default  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$)
$(info $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$)
$(info $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$)
$(info $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$)
endif

ifeq ($(LD_FLAGS),)
$(info LD_FLAGS not set)
endif

#BLAS_PATH=/home/abduljm/libs/petsc-3.9.3/arch-linux2-c-debug
#BLAS_LIB=${BLAS_PATH}/lib
#
#LAPACK_PATH=/home/abduljm/libs/petsc-3.9.3/arch-linux2-c-debug
#LAPACK_LIB=${LAPACK_PATH}/lib
#
#LIB_BLAS=-L${BLAS_LIB} -lblas
#LIB_LAPACK=-L${LAPACK_LIB} -llapack 
#
#LD_FLAGS+=${LIB_LAPACK} ${LIB_BLAS}

ifneq (,$(findstring intel,$(ACOUS_ENV)))
F77=fn
F90=fn
CXX=mpiicpc
##CXX=CC
ifeq ($(debug),1)
                F77FLAGS=-g -CB
                F90FLAGS=-g -CB
                CXXFLAGS+=-g
else
                F77FLAGS=-O3
                F90FLAGS=-O3
                CXXFLAGS+=-O3 -fabi-version=6 -funroll-loops -qopenmp -mavx 
endif   
		LD_FLAGS+=-qopenmp 

else ifneq (,$(findstring cray,$(ACOUS_ENV)))
F77=ftn
F90=ftn
CXX=CC
ifeq ($(debug),1)
                F77FLAGS=-g -CB
                F90FLAGS=-g -CB
                CXXFLAGS+=-g
else
                F77FLAGS=-O3
                F90FLAGS=-O3
                #CXXFLAGS+=-O3 -qopenmp -axCORE-avx2
                CXXFLAGS+=-O3 -fopenmp -mavx
endif   
LD_FLAGS+=-fopenmp
#LD_FLAGS+=-L/opt/intel/composer_xe_2015.2.164/tbb/lib/intel64/gcc4.1/ -Wl,-R/opt/intel/composer_xe_2015.2.164/tbb/lib/intel64/gcc4.1/ -ltbb -dynamic
#LD_FLAGS+=-L/opt/intel/composer_xe_2015.2.164/compiler/lib/intel64/ -Wl,-R/opt/intel/composer_xe_2015.2.164/compiler/lib/intel64/ -lifcore -lsvml -lirng -lcilkrts -lirc -dynamic

else 
F77=gfortran
F90=gfortran
CXX=mpicxx
ifeq ($(debug),1)
                F77FLAGS=-g -fbounds-check -fbacktrace # -g #
                F90FLAGS=-g -Wunused-variable -Wunused-parameter -Wunused-dummy-argument -fbounds-check -fbacktrace # -g #
                CXXFLAGS+=-g -fabi-version=6 -Wunused-variable -funroll-loops -fopenmp 
else
                F77FLAGS=-O3 -fbounds-check -fbacktrace # -g #
                F90FLAGS=-O3 -fbounds-check -fbacktrace # -g #
                CXXFLAGS+=-O3 -fabi-version=6 -ffast-math -funroll-loops -fopenmp -mavx
endif
		LD_FLAGS+=-fopenmp -lgfortran -lpetsc 
endif



#LIB_EXAFMM_ACOUSTICS=${EXAFMM_PATH}/wrappers/libacoustics_exafmm_helmholtz.a
ACOUSTICS_INCLUDE=-Iinclude/ -I${PETSC_DIR}/include/ -I${PETSC_DIR}/${PETSC_ARCH}/include/ -Iinclude/fmm -I${PARMETIS_PATH}/metis/include -I${PARMETIS_PATH}/include
ifneq (,$(findstring tbb,$(THREADING_MODEL)))
ACOUSTICS_INCLUDE+=-I${TBBROOT}/include/ 
endif
LIBS=${LIB_EXAFMM_ACOUSTICS} ${PETSC_LIB} ${LIB_PARMETIS} ${LD_FLAGS}
 
F90SRC= fsrc/global_com.f90 fsrc/gaussintegraldata.f90 fsrc/global_geom.f90 \
				fsrc/global_dim.f90  fsrc/curve_para.f90 fsrc/inmom.f90 fsrc/incurve.f90\
				fsrc/soft_nys.f90 fsrc/xyz_position.f90 fsrc/inci_rhs.f90 fsrc/gmres_fannew.f90\
				fsrc/rcs.f90 fsrc/near_scafield.f90 fsrc/solver.f90

F90OBJ= fsrc/global_com.o fsrc/gaussintegraldata.o fsrc/global_geom.o \
				fsrc/global_dim.o  fsrc/curve_para.o fsrc/inmom.o fsrc/incurve.o\
				fsrc/soft_nys.o fsrc/xyz_position.o fsrc/inci_rhs.o fsrc/gmres_fannew.o\
				fsrc/rcs.o fsrc/near_scafield.o fsrc/solver.o				

CPP=bemfmm_test_mpi.cxx
FTN=fbin/frequency.f90
MOT=$(F90SRC:.f90=.o)
EXE=$(CPP:.cxx=.o)
FEXE=$(FTN:.cxx=.o)


.SUFFIXES:      .o .f .c .C .h .f90 .cxx

.f90.o:
	$(F90) $(F90FLAGS) -c -o $@ $<
.cxx.o:
	$(CXX) $(CXXFLAGS) $(ACOUSTICS_INCLUDE) -c $<

.DEFAULT_GOAL := bemfmm_test_mpi

#exafmm_lib:
#	cd ${EXAFMM_PATH}/wrappers/ && $(MAKE) clean && $(MAKE) libacoustics_exafmm_helmholtz.a -j

bemfmm_test_mpi: $(EXE) 
	$(CXX) $(EXE) $(LIBS) $(LINK_FORT_ACOUSTICS) -o $@ 

	
cleanbin:
	rm -f bemfmm_test_mpi
	rm -f *.o
	rm -f *.out
	
bin: cleanbin bemfmm_test_mpi

clean::
	rm -f *.o
	rm -f fsrc/*.o
	rm -f *.out
	rm -f *.dat
	rm -f fort.*
	rm -f *.time
	rm -f *.mod
	rm -f *.o.optdbg
	rm -f *.bak
	rm -f *~	
	rm -f bemfmm_test_mpi
