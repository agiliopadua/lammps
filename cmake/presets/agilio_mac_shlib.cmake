# preset that turns on a wide range of packages, some of which require
# external libraries. Compared to all_on.cmake some more unusual packages
# are removed. The resulting binary should be able to run most inputs.

set(ALL_PACKAGES
  COLVARS
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
  PYTHON
  RIGID
  USER-MISC
  VORONOI)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(BUILD_TOOLS ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)
set(CMAKE_INSTALL_PREFIX "/Users/apadua" CACHE STRING "" FORCE)
set(LDFLAGS="-L/opt/homebrew/opt/libomp/lib" CACHE STRING "" FORCE)
set(CPPFLAGS="-I/opt/homebrew/opt/libomp/include" CACHE STRING "" FORCE)
