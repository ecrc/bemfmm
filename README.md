# README # [![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

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


