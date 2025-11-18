cmake_minimum_required(VERSION 3.19...3.30)

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

file(GENERATE OUTPUT .gitignore CONTENT "*")

set(CMAKE_CXX_STANDARD 17)

# --- user options

include(cmake/options.cmake)
message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}  install prefix: ${CMAKE_INSTALL_PREFIX}  CMake ${CMAKE_VERSION}  Toolchain ${CMAKE_TOOLCHAIN_FILE}")

# --- external libs

if(FCLAW_ENABLE_CUDACLAW)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  find_package(CUDAToolkit REQUIRED)
endif()

# --- p4est, sc
include(cmake/GitSubmodule.cmake)

if(FCLAW_ENABLE_SUBMODULES)

  git_submodule(${CMAKE_CURRENT_SOURCE_DIR}/p4est)
  set(P4EST_BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(SC_BUILD_TESTING OFF CACHE BOOL "" FORCE)
  add_subdirectory(p4est)

else()

  find_package(P4EST CONFIG REQUIRED)

endif()

# --- ThunderEgg
if(FCLAW_ENABLE_THUNDEREGG)
  if(FCLAW_ENABLE_SUBMODULES)

    git_submodule(${CMAKE_CURRENT_SOURCE_DIR}/thunderegg)
    add_subdirectory(thunderegg)

  else()

    find_package(ThunderEgg REQUIRED)

  endif()
endif()

# --- configure

include(cmake/config.cmake)
include(cmake/compilers.cmake)

# Conda's HDF5 never works with the system compiler
if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_IGNORE_PREFIX_PATH $ENV{CONDA_PREFIX})
  list(APPEND CMAKE_IGNORE_PATH $ENV{CONDA_PREFIX}/bin)
  # need CMAKE_IGNORE_PATH for CMake < 3.23
  # and to ensure system env var PATH doesn't interfere despite CMAKE_IGNORE_PREFIX_PATH
endif()

if(FCLAW_ENABLE_HDF5 AND NOT TARGET HDF5::HDF5)
  if(FCLAW_ENABLE_SUBMODULES)

    include(cmake/hdf5.cmake)

  else()

    find_package(HDF5 REQUIRED)

  endif()
endif()



# --- test
if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)


# --- ForestClaw
add_subdirectory(src)

# --- Clawpatch
add_subdirectory(src/patches/clawpatch)


# --- Clawpack 4.6
if(FCLAW_ENABLE_CLAWPACK)

  add_subdirectory(src/solvers/fc2d_clawpack4.6)
  add_subdirectory(src/solvers/fc2d_clawpack5)
  add_subdirectory(src/solvers/fc3d_clawpack46)

endif(FCLAW_ENABLE_CLAWPACK)



# --- geoclaw
if(FCLAW_ENABLE_GEOCLAW)

  add_subdirectory(src/solvers/fc2d_geoclaw)

endif(FCLAW_ENABLE_GEOCLAW)



# --- cudaclaw
if(FCLAW_ENABLE_CUDACLAW)

  add_subdirectory(src/solvers/fc2d_cudaclaw)

endif(FCLAW_ENABLE_CUDACLAW)



# --- ThunderEgg
if(FCLAW_ENABLE_THUNDEREGG)

  add_subdirectory(src/solvers/fc2d_thunderegg)
  add_subdirectory(src/solvers/fc3d_thunderegg)

endif(FCLAW_ENABLE_THUNDEREGG)
# --- package


# -- put together in final library

include(cmake/pkgconf.cmake)
include(cmake/install.cmake)


# --- examples
if(FCLAW_ENABLE_APPLICATIONS)
  add_subdirectory(applications)
endif()

include(cmake/documentation.cmake)
