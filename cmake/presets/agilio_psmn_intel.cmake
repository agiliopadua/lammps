# preset that turns on a wide range of packages, some of which require
# external libraries. Compared to all_on.cmake some more unusual packages
# are removed. The resulting binary should be able to run most inputs.

set(ALL_PACKAGES
  CORESHELL
  DRUDE
  ELECTRODE
  FEP
  KSPACE
  MANYBODY
  MC
  MISC
  MOFFF
  MOLECULE
  OPENMP
  PLUMED
  RIGID
  USER-MISC
  VORONOI)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(BUILD_TOOLS ON CACHE BOOL "" FORCE)
set(CMAKE_INSTALL_PREFIX "/home/apadua" CACHE STRING "" FORCE)

set(CMAKE_CXX_COMPILER "icpc" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "icc" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "" FORCE)
set(MPI_CXX "icpc" CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
unset(HAVE_OMP_H_INCLUDE CACHE)
set(OpenMP_C "icc" CACHE STRING "" FORCE)
set(OpenMP_C_FLAGS "-qopenmp -qopenmp-simd" CACHE STRING "" FORCE)
set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "" FORCE)
set(OpenMP_CXX "icpc" CACHE STRING "" FORCE)
set(OpenMP_CXX_FLAGS "-qopenmp -qopenmp-simd" CACHE STRING "" FORCE)
set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
set(OpenMP_Fortran_FLAGS "-qopenmp -qopenmp-simd" CACHE STRING "" FORCE)
set(OpenMP_omp_LIBRARY "libiomp5.so" CACHE PATH "" FORCE)
