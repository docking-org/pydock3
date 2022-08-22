# GNU gfortran linux makefile header

# Work in progres...

ifeq ($(COMPILER), gnu)
F77 = gfortran
CC = gcc
GPROF = gprof

F_VERSION=f2003

# === OPTIMIZATION FLAGS ===
# Optimizations for specific 32b/64b targets
OPT64 = 
OPT32 = 
ifndef DEBUG
# -- FULL optimizations
OPT =  -O3

# -- Optimizations for profiling
# -- Turn off inlining so function boundaries are visible, though this will reduce performance
OPTPROF =  -O2
else
# -- DEBUG optimizations
OPT = -C -fbacktrace -fbounds-check -Wall -g 
#OPT = -C -fbacktrace -g 

# -- Optimizations for profiling
# -- Turn off inlining so function boundaries are visible, though this will reduce performance
OPTPROF =  -O0 
endif
# ==========================


 
# === HW / ARCHITECTURE OPTIONS ===
# Set ARCH32 and ARCH64 to the highest architecture setting possible that 
# is the lowest common denominator of all the 64b (or 32b) servers in the 
# server grid.
ARCH64 = -m64
ARCH32 = -m32
# =================================


# === COMPILER FLAGS ===
#FFLAGS_COMMON = -fconvert=swap -fcray-pointer -std=f95 -fmax-errors=2 -fall-intrinsics  $(ARCH) 
#FFLAGS_COMMON =  -std=f95 -fcray-pointer  $(ARCH) 
# Intrnsics needed for built-ins
FFLAGS_COMMON = -std=$(F_VERSION) $(ARCH) -fall-intrinsics -fconvert=swap
CFLAGS_COMMON = $(ARCH)  
FFLAGS =  $(FFLAGS_COMMON) -g $(OPT) -static
CFLAGS =  $(CFLAGS_COMMON) -g $(OPT)
# FFLAGS for Profiling.  Add in -g symbols
FFLAGS_PROF =  $(FFLAGS_COMMON) $(OPTPROF) -g -Mprof=func
# GNU gprof
FFLAGS_GPROF = $(FFLAGS_COMMON) $(OPTPROF) -pg  
# Intel VTune
FFLAGS_VTUNE =  $(FFLAGS_COMMON) $(OPTPROF) -g -Mdwarf2
CFLAGS_VTUNE =  $(CFLAGS_COMMON) $(OPTPROF) -g -Mdwarf2 
# ======================

# === LINKER OPTIONS ===
# Linker libraries
FFLIBS = -lfgz$(SIZE) -lz
#-static-libgfortran
# ======================
endif
