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
  MOLECULE
  OPENMP
  RIGID
  USER-MISC)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(BUILD_TOOLS ON CACHE BOOL "" FORCE)
set(CMAKE_CXX_COMPILER "g++" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "gcc" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "" FORCE)
set(CMAKE_INSTALL_PREFIX "/home/apadua" CACHE STRING "" FORCE)

