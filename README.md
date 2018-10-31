# BEMFMM #

An extreme-scale Fast Multipole Method (FMM)-accelerated Boundary Integral Equation (BIE) solver for wave scattering. The solver name stands for Boundary Element Method (BEM) combined with FMM (BEMFMM). The application code calculates scattered field due to an excitation source at a specific point in space or infinity. The Krylov subspace liner solver is featured by GMRES iterative method, inside which FMM is used to implement the Matrix-Vector multiplication (MatVec) kernel. The solver is highly optimized for both shared- and distributed-memory architectures, and support optimal architecture-specific and algorithm-aware partitioning, load balancing, and communication reducing mechanisms. The solver code utilizes both *ExaFMM* and *FMMLIB3D*. Followings are two diagrams that depict the underlying implementation of the solver code.

![Image of BEMFMM workflow](https://github.com/ecrc/BEMFMM/blob/master/img/workflow.png)

![Image of the implemented FMM](https://github.com/ecrc/BEMFMM/blob/master/img/fmmG.png)

### What is this repository for? ###

This repository is a large-scale FMM-based Boundary Element Solver that calculates scattered field due to an exciation source at a specific point in space or infinity. 


### Version ###
1.0

### Prerequisites ###
You need at least a Fortran 90 Compiler, libraries needed for now are as follows:

- MPI >= 2.0 e.g. OpenMPI, MPICH

- Lapack and Blas

- Intel TBB

- ParMetis


## Running a test case ###
After you have installed the required libraries, do

make clean

make 


