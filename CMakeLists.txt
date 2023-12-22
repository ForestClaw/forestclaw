cmake_minimum_required(VERSION 3.15...3.28)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "please use out-of-source build
  cmake -Bbuild
  ")
endif()

include(cmake/GitVersionGen.cmake)

project(FORESTCLAW
  DESCRIPTION "Quadtree/octree adaptive PDE solver based based on p4est."
  HOMEPAGE_URL https://github.com/ForestClaw/forestclaw
  LANGUAGES C CXX Fortran
  VERSION ${FORESTCLAW_VERSION}
)

include(CTest)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

set(CMAKE_CXX_STANDARD 17)

# --- user options

include(cmake/options.cmake)
message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}  install prefix: ${CMAKE_INSTALL_PREFIX}  CMake ${CMAKE_VERSION}  Toolchain ${CMAKE_TOOLCHAIN_FILE}")

# --- external libs

if(mpi)
  find_package(MPI COMPONENTS C CXX REQUIRED)
endif()

if(cudaclaw)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  find_package(CUDAToolkit REQUIRED)
endif()

#find_package(LAPACK)
#find_package(BLAS)
find_package(ZLIB)

# --- p4est, sc

include(cmake/sc.cmake)
include(cmake/p4est.cmake)

# --- ThunderEgg
if(thunderegg)
  include(cmake/ThunderEgg.cmake)
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

include(cmake/documentation.cmake)
