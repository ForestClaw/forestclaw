option(mpi "use MPI library")
option(openmp "use OpenMP")
option(applications "build applications" ON)

option(clawpack "build Clawpack")
option(geoclaw "build Geoclaw")
option(cudaclaw "build CudaClaw")
option(thunderegg "build ThunderEgg")

option(thunderegg_external "force build of ThunderEgg")
option(p4est_external "force build of p4est")
option(sc_external "force build of libsc")

set(CMAKE_TLS_VERIFY on)

# --- default install directory under build/local
# users can specify like "cmake -B build -DCMAKE_INSTALL_PREFIX=~/mydir"
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # will not take effect without FORCE
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/local" CACHE PATH "Install top-level directory" FORCE)
endif()

# enable needed dependencies
if(clawpack)
  set(clawpatch ON)
  set(clawpack4.6 ON)
  set(clawpack5 ON)
  set(clawpatch ON)
  set(clawpack3_46 ON)
endif(clawpack)

if(geoclaw)
  set(clawpatch ON)
endif(geoclaw)

if(cudaclaw)
  set(clawpatch ON)
  set(clawpack4.6 ON)
  set(clawpack5 ON)
endif(cudaclaw)

if(thunderegg)
  set(clawpatch ON)
  set(clawpack4.6 ON)
endif(thunderegg)

# Rpath options necessary for shared library install to work correctly in user projects
set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

# Necessary for shared library with Visual Studio / Windows oneAPI
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS true)

# --- auto-ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()
