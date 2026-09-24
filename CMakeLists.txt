# Builds the programs blastermaster runs (sources in native/) and installs them into the
# Python package, at pydock3/blastermaster/bin/.
cmake_minimum_required(VERSION 3.15)
project(pydock3_programs LANGUAGES C CXX Fortran)

set(BIN_DIR pydock3/blastermaster/bin)
option(PYDOCK3_STATIC_RUNTIME "Link the compiler runtime libraries statically (for wheels)" OFF)

# Same results on every platform: no fused multiply-adds (used on ARM CPUs), and no
# vectorized Fortran math functions (glibc's libmvec, less accurate than the scalar ones).
# The Fortran is FORTRAN 77 that, like the PGI compilers it was written for, expects local
# variables to be static (so zero-initialized and kept between calls).
add_compile_options(
  -ffp-contract=off
  "$<$<COMPILE_LANGUAGE:Fortran>:-std=legacy;-fno-automatic;-fno-tree-vectorize>"
)

if(PYDOCK3_STATIC_RUNTIME)
  if(WIN32)
    add_link_options(-static)
  else()
    add_link_options(
      -static-libgcc
      $<$<LINK_LANGUAGE:CXX>:-static-libstdc++>
      $<$<LINK_LANGUAGE:Fortran>:-static-libgfortran>
    )
  endif()
endif()

# --- Fortran programs from DOCK (originally built with PGI compilers) ---------------------

add_executable(filt native/filt/filter1.f)
target_compile_options(filt PRIVATE -ffixed-line-length-none)

add_executable(sphgen native/sphgen/sphgen.f)

add_executable(chemgrid
  native/chemgrid/chemgrid.f native/chemgrid/dconst.f native/chemgrid/ddist.f
  native/chemgrid/dist.f native/chemgrid/grdout.f native/chemgrid/parmrec.f
)
# DOCK reads the grid as big-endian; the atom-type hash relies on integer overflow wrapping
target_compile_options(chemgrid PRIVATE -fconvert=big-endian -fwrapv)

add_executable(solvmap
  native/solvmap/solvmap.f native/solvmap/sas.f native/solvmap/cube.f native/solvmap/memalloc.c
)
target_compile_options(solvmap PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fcray-pointer>)

file(GLOB QNIFFT_SOURCES CONFIGURE_DEPENDS native/qnifft/*.f)
add_executable(qnifft ${QNIFFT_SOURCES})
# DOCK reads the phi map as big-endian
target_compile_options(qnifft PRIVATE -ffixed-line-length-132 -fconvert=big-endian)

# --- dms (UCSF Computer Graphics Lab): client and its compute server ----------------------

add_executable(dms
  native/dms/compute.c native/dms/fwritev.c native/dms/input.c native/dms/ms.c
  native/dms/output.c native/dms/emalloc.c native/dms/tokenize.c native/dms/dms_paths.c
  native/dms/libpdb/pdb_read.c native/dms/libpdb/pdb_sprntf.c native/dms/libpdb/pdb_sscanf.c
  native/dms/libpdb/pdb_write.c native/dms/libpdb/pdbrun.c native/dms/libpdb/ms.c
)
target_include_directories(dms PRIVATE native/dms native/dms/libpdb)

add_executable(dmsd
  native/dms/dmsd/server.c native/dms/dmsd/viewat.c native/dms/dmsd/lookat.c
  native/dms/dmsd/freadv.c native/dms/dmsd/emalloc.c
)
target_include_directories(dmsd PRIVATE native/dms)
# dmsd uses correctly rounded acos, atan2, sin, cos (CORE-MATH), so that dms computes the same
# surface on every platform (the platforms' math libraries can differ in the last bit)
add_library(core_math STATIC
  native/core-math/acos.c native/core-math/atan2.c native/core-math/sin.c native/core-math/cos.c
)
target_link_libraries(dmsd PRIVATE core_math)
if(NOT WIN32)
  target_link_libraries(dmsd PRIVATE m)
endif()

# K&R C: implicit declarations are valid C89
foreach(target dms dmsd)
  target_compile_options(${target} PRIVATE -std=gnu89 -w)
endforeach()

# --- reduce (Richardson Lab), without its Python bindings --------------------------------

set(CMAKE_CXX_STANDARD 11)
add_subdirectory(native/reduce/libpdb)
add_subdirectory(native/reduce/toolclasses)
add_executable(reduce
  native/reduce/reduce_src/main.cpp native/reduce/reduce_src/GraphToHoldScores.cpp
  native/reduce/reduce_src/reduce.cpp native/reduce/reduce_src/CTab.cpp
  native/reduce/reduce_src/ElementInfo.cpp native/reduce/reduce_src/StdResH.cpp
  native/reduce/reduce_src/ResBlk.cpp native/reduce/reduce_src/AtomConn.cpp
  native/reduce/reduce_src/AtomPositions.cpp native/reduce/reduce_src/DotSph.cpp
  native/reduce/reduce_src/Mover.cpp native/reduce/reduce_src/Rot.cpp
  native/reduce/reduce_src/Rot3Fold.cpp native/reduce/reduce_src/RotAromMethyl.cpp
  native/reduce/reduce_src/RotDonor.cpp native/reduce/reduce_src/FlipMemo.cpp
  native/reduce/reduce_src/CliqueList.cpp native/reduce/reduce_src/AtomDescr.cpp
  native/reduce/reduce_src/PDBrec.cpp native/reduce/reduce_src/MoveableNode.cpp
  native/reduce/reduce_src/hybrid_36_c.c
)
target_compile_definitions(reduce PRIVATE
  AROMATICS_ACCEPT_HBONDS LEFT_JUSTIFY_NUC_RES_OK
  HET_DICTIONARY="/usr/local/reduce_wwPDB_het_dict.txt"  # unused: pydock3 passes -db
  HET_DICTOLD="/usr/local/reduce_het_dict.txt"
)
target_include_directories(reduce PRIVATE native/reduce/libpdb native/reduce/toolclasses)
target_link_libraries(reduce PRIVATE pdb++ toolclasses)

install(TARGETS filt sphgen chemgrid solvmap qnifft dms dmsd reduce RUNTIME DESTINATION ${BIN_DIR})
