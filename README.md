# README #

### What is this repository for? ###

The repository has a Fortran code that solves the 3D Helmholtz equation to find the acoustic wave scattering by near and far field objects on a target surface
Nystrom method is used for discretizing the resulting Boundary Integral Equations. 
Eventually FMM will be used to do the inner Mat-Vec products for far field interactions.


### Version ###
1.0

### Prerequisites ###
You need at least a Fortran 90 Compiler, libraries needed for now are as follows:

- An MPI 2.0 library e.g. OpenMPI, MPICH

- Lapack

- BLAS

- FFTW3

## Running a test case ###
After you have installed the required libraries, do

make clean

make 

### Note that this repository will be gradually changed from Fortran to C++ for future reusability and potential performance enhancements ###

