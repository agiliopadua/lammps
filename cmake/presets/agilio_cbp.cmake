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
set(CMAKE_INSTALL_PREFIX "/projects/DepartementChimie/lammps" CACHE STRING "" FORCE)
set(DOWNLOAD_PLUMED "no" CACHE STRING "" FORCE)
set(PLUMED_MODE "shared" CACHE STRING "" FORCE)
set(PKG_CONFIG_PATH "/projects/DepartementChimie/plumed/lib/pkgconfig" CACHE STRING "" FORCE)
