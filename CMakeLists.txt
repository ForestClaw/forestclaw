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

include(FetchContent)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

file(GENERATE OUTPUT .gitignore CONTENT "*")

set(CMAKE_CXX_STANDARD 17)

# --- user options

include(cmake/options.cmake)
message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}  install prefix: ${CMAKE_INSTALL_PREFIX}  CMake ${CMAKE_VERSION}  Toolchain ${CMAKE_TOOLCHAIN_FILE}")

# --- external libs

if(cudaclaw)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  find_package(CUDAToolkit REQUIRED)
endif()

# --- p4est, sc
include(cmake/GitSubmodule.cmake)

git_submodule(${CMAKE_CURRENT_SOURCE_DIR}/sc)
FetchContent_Declare(SC SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sc)
FetchContent_MakeAvailable(SC)

git_submodule(${CMAKE_CURRENT_SOURCE_DIR}/p4est)
FetchContent_Declare(P4EST SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/p4est)
FetchContent_MakeAvailable(P4EST)

# --- ThunderEgg
if(thunderegg)
  include(cmake/ThunderEgg.cmake)
endif(thunderegg)

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
if(fclaw_applications)
  add_subdirectory(applications)
endif()

include(cmake/documentation.cmake)
