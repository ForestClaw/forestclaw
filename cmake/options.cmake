option(mpi "use MPI library")
option(openmp "use OpenMP")
option(applications "build applications" ON)

option(clawpatch "build Clawpatch")
option(clawpack "build Clawpack")
option(geoclaw "build Geoclaw")
option(cudaclaw "build CudaClaw")
option(thunderegg "build ThunderEgg")

option(thunderegg_external "force build of ThunderEgg")
option(p4est_external "force build of p4est")
option(sc_external "force build of libsc")

option(CMAKE_TLS_VERIFY "verify HTTPS certs" on)

# --- default install directory under build/local
# users can specify like "cmake -B build -DCMAKE_INSTALL_PREFIX=~/mydir"
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND PROJECT_IS_TOP_LEVEL)
  # will not take effect without FORCE
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/local" CACHE PATH "Install top-level directory" FORCE)
endif()

unset(${PROJECT_NAME}_stdfs_link_flags)
if( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.1.0") OR
  (LINUX AND CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "23") OR
  (CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "23.11") )
set(${PROJECT_NAME}_stdfs_link_flags stdc++fs stdc++)
endif()
# GCC < 9.1 needs -lstdc++ to avoid C main program link error
# NVHPC at least 23.11 and newer doesn't need the flags, but at least 23.5 and older do.
# INtel oneAPI 2021.1 and older needs, but 2023 and newer doesn't. (not sure about 2022)


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

# Necessary for shared library with Visual Studio / Windows oneAPI
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS true)
