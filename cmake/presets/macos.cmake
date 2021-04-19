# preset that turns on a wide range of packages, some of which require
# external libraries. Compared to all_on.cmake some more unusual packages
# are removed. The resulting binary should be able to run most inputs.

set(ALL_PACKAGES PYTHON CLASS2 CORESHELL DIPOLE KSPACE MANYBODY MC MISC
        MOLECULE QEQ RIGID VORONOI
        USER-COLVARS USER-DRUDE USER-FEP USER-MISC USER-MOFFF USER-MOLFILE
        USER-OMP USER-REACTION)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

set(BUILD_TOOLS ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)
set(CMAKE_INSTALL_PREFIX "/Users/apadua" CACHE STRING "" FORCE)
