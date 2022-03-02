# preset that turns on a wide range of packages, some of which require
# external libraries. Compared to all_on.cmake some more unusual packages
# are removed. The resulting binary should be able to run most inputs.

set(ALL_PACKAGES
  COLVARS
  CORESHELL
  DRUDE
  FEP
  KSPACE
  MC
  MISC
  MOFFF
  MOLECULE
  MOLFILE
  OPENMP
  RIGID
  USER-MISC
  VORONOI)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(CMAKE_INSTALL_PREFIX "/home/apadua" CACHE STRING "" FORCE)
