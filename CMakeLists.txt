cmake_minimum_required(VERSION 3.13...3.21)

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
find_package(ZLIB REQUIRED)

# --- find p4est or build if not found

if(NOT sc_external)
  if(NOT SC_ROOT)
    if(THUNDEREGG_ROOT)
      set(SC_ROOT ${THUNDEREGG_ROOT})
    elseif(P4EST_ROOT)
      set(SC_ROOT P4EST_ROOT)
    endif()
  endif()
  find_package(SC)
endif()
if(NOT TARGET SC::SC)
  include(cmake/sc.cmake)
endif()

# --- find libsc or build if not found

if(NOT p4est_external)
  if(NOT P4EST_ROOT AND DEFINED THUNDEREGG_ROOT)
    set(P4EST_ROOT ${THUNDEREGG_ROOT})
  endif()
  find_package(P4EST)
endif()
if(NOT TARGET P4EST::P4EST)
  include(cmake/p4est.cmake)
endif()

# --- ThunderEgg
if(thunderegg)
  if(NOT thunderegg_external)
    if(NOT ThunderEgg_ROOT AND DEFINED THUNDEREGG_ROOT)
      set(ThunderEgg_ROOT ${THUNDEREGG_ROOT})
    endif()
    find_package(ThunderEgg)
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
if(clawpack3)

  add_subdirectory(src/solvers/fc3d_clawpack46)

endif(clawpack3)


include(cmake/pkgconf.cmake)
include(cmake/install.cmake)

# --- examples
add_subdirectory(applications)
