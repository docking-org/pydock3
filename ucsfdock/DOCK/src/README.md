Compiling DOCK
==============

**Note**: all commands given from this (DOCK/docking/DOCK/src) directory.

Environmental Variables
-----------------------
 - **COMPILER**: Which compiler suite to use. Options: pgf (default), gcc, ifort
 - **SIZE**: Which architecture to generate. Options: 64 (default), 32

Preparation
-----------
Before compiling ensure that you have removed all intermediate files from other compilers.

    cd libfgz
    make clean
    cd ../i386
    make clean
    cd ..

Compiling With PGI (Default)
----------------------------
The PGI compilation of DOCK is heavily optimized and will result in the best performance by far. This is the recommended 
method of compilation and how the official binaries are generated.

This method requires a functional installation of the Portland Group's Fortran and C compilers!

### Dependencies (CentOS)
 - PGI Workstation (http://www.pgroup.com/products/pgiworkstation.htm)
  * PGFORTRAN
  * PGCC

### Compilation

    cd libfgz
    COMPILER=pgf SIZE=64 make
    cd ../i386
    COMPILER=pgf SIZE=64 make


Compiling With GNU
------------------
The GNU compilation method is provided for completeness and flexibility. No guarantees of speed are made.

### Dependencies (CentOS)
 - gcc
 - gcc-gfortran
 - glibc-devel 
 - glibc-static
 - zlib-devel 
 - zlib-static
  
### Compilation

    cd libfgz
    COMPILER=gnu SIZE=64 make
    cd ../i386
    COMPILER=gnu SIZE=64 make


Compiling With Intel
--------------------
**WARNING:** This method is incomplete and has not been tested quite a while!

### Dependencies (CentOS):
 - Intel Compiler Collection
  * IFort

### Compiliation

    cd libfgz
    COMPILER=ifort SIZE=64 make
    cd ../i386
    COMPILER=ifort SIZE=64 make
