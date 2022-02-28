cmake_minimum_required(VERSION 3.13...3.20)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 70)
endif()

project(FORESTCLAW
  DESCRIPTION "Quadtree/octree adaptive PDE solver based based on p4est."
  HOMEPAGE_URL https://github.com/ForestClaw/forestclaw
  LANGUAGES C CXX Fortran
  VERSION 0.1.4880)

include(CTest)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

set(CMAKE_CXX_STANDARD 17)
# --- user options

include(cmake/options.cmake)
message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}  install prefix: ${CMAKE_INSTALL_PREFIX}")

# --- external libs

if(mpi)
  find_package(MPI COMPONENTS C CXX REQUIRED)
endif()

if(cudaclaw)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
endif()

#find_package(LAPACK)
#find_package(BLAS)
find_package(ZLIB)

# --- p4est, sc

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

if(NOT sc_external)
  find_package(SC HINTS ${THUNDEREGG_ROOT} ${P4EST_ROOT})
endif()
if(NOT TARGET SC::SC)
  include(cmake/sc.cmake)
endif()


if(NOT p4est_external)
  find_package(P4EST HINTS ${THUNDEREGG_ROOT})
endif()
if(NOT TARGET P4EST::P4EST)
  include(cmake/p4est.cmake)
endif()



# --- ThunderEgg
if(thunderegg)
  if(NOT thunderegg_external)
    find_package(ThunderEgg HINTS ${THUNDEREGG_ROOT})
  endif()
  if(NOT TARGET ThunderEgg::ThunderEgg)
    include(cmake/ThunderEgg.cmake)
  endif()
endif(thunderegg)


if(cudaclaw)
endif(cudaclaw)
# --- configure


include(cmake/config.cmake)
include(cmake/compilers.cmake)

# --- test
if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)


# --- ForestClaw
add_subdirectory(src)



# --- Clawpatch
if(clawpatch)

  add_subdirectory(src/patches/clawpatch)
  
endif(clawpatch)



# --- Clawpack 4.6
if(clawpack4.6)

  add_subdirectory(src/solvers/fc2d_clawpack4.6)

endif(clawpack4.6)



# --- Clawpack 5
if(clawpack5)

  add_subdirectory(src/solvers/fc2d_clawpack5)

endif(clawpack5)



# --- geoclaw
if(geoclaw)

  add_subdirectory(src/solvers/fc2d_geoclaw)

endif(geoclaw)



# --- cudaclaw
if(cudaclaw)

  add_subdirectory(src/solvers/fc2d_cudaclaw)

endif(cudaclaw)



# --- ThunderEgg
if(thunderegg)

  add_subdirectory(src/solvers/fc2d_thunderegg)

endif(thunderegg)
# --- package


# --- clawpack3
if(clawpack3_46)

  add_subdirectory(src/solvers/fc3d_clawpack46)

endif(clawpack3_46)


include(cmake/pkgconf.cmake)
include(cmake/install.cmake)


# --- examples
if(applications)
  add_subdirectory(applications)
endif(applications)


