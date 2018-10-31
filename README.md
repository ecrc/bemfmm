# BEMFMM #

An extreme-scale Fast Multipole Method (FMM)-accelerated Boundary Integral Equation (BIE) solver for wave scattering. The solver name stands for Boundary Element Method (BEM) combined with FMM (BEMFMM). The application code calculates scattered field due to an excitation source at a specific point in space or infinity. The Krylov subspace liner solver is featured by GMRES iterative method, inside which FMM is used to implement the Matrix-Vector multiplication (MatVec) kernel. The solver is highly optimized for both shared- and distributed-memory architectures, and support optimal architecture-specific and algorithm-aware partitioning, load balancing, and communication reducing mechanisms. The solver code utilizes two state-of-the-art FMM implementation for oscillatory kernels, namely *ExaFMM* and *FMMLIB3D*. Followings are two diagrams that depict the underlying implementation of the solver code.

![Image of BEMFMM workflow](https://github.com/ecrc/BEMFMM/blob/master/img/workflow.png = 128x128)

![Image of the implemented FMM](https://github.com/ecrc/BEMFMM/blob/master/img/fmmG.png = 128x128)

### Requirements ###

* C/C++ Compiler (e.g., GNU Compiler -- https://www.gnu.org/software/gcc/)
* MPI  (e.g., MPICH -- http://www.mpich.org/)
* LAPACK (e.g., NETLIB LAPACK -- http://www.netlib.org/lapack/)
* Intel TBB (https://www.threadingbuildingblocks.org/)
* ParMETIS (http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)

The repository includes LAPACK, TBB, and ParMETIS. Therefore, you may not need to install yours, use can just use the ones that are included herein. However, if you have a better implementation that you wish to link to, then you can just install it on your software environment, and directly link to your existing implementation. Hence, the minimum requirements to run BEMFMM are:

* C/C++ Compiler (e.g., GNU Compiler -- https://www.gnu.org/software/gcc/)
* MPI  (e.g., MPICH -- http://www.mpich.org/)

Please have these two dependencies configured and installed on your system before running the solver code.

## Running a test case ###
After you have installed the required libraries, do

make clean

make 


