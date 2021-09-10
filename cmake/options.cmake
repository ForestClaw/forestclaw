option(mpi "use MPI library")
option(openmp "use OpenMP")

option(clawpack "build Clawpack")
option(geoclaw "build Geoclaw")
option(cudaclaw "build CudaClaw")
option(thunderegg "build ThunderEgg")
option(clawpack3 "build Clawpack 3")

option(thunderegg_external "force build of ThunderEgg")
option(p4est_external "force build of p4est")
option(sc_external "force build of libsc")

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

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
  set(clawpack3_46 ON)
endif(clawpack)

if(clawpack3)
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