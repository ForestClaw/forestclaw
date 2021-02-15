cmake_minimum_required(VERSION 3.13)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Debug or Release")
endif()

project(ForestClawApps LANGUAGES C CXX Fortran)

include(CTest)

find_package(MPI COMPONENTS C Fortran REQUIRED)
find_package(ZLIB REQUIRED)

# --- get installed ForestClaw package
find_package(FORESTCLAW CONFIG REQUIRED)

if(EXISTS ${FORESTCLAW_DIR}/../SC/SCTargets.cmake)
  include(${FORESTCLAW_DIR}/../SC/SCTargets.cmake)

  add_library(SC::SC INTERFACE IMPORTED)
  target_link_libraries(SC::SC INTERFACE SC::sc)
endif()

include(${FORESTCLAW_DIR}/../P4EST/P4ESTTargets.cmake)
add_library(P4EST::P4EST INTERFACE IMPORTED)
target_link_libraries(P4EST::P4EST INTERFACE P4EST::p4est)

include(${FORESTCLAW_DIR}/FORESTCLAWTargets.cmake)
add_library(FORESTCLAW::FORESTCLAW INTERFACE IMPORTED)
target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE FORESTCLAW::forestclaw)

# --- apps

add_subdirectory(clawpack/advection/2d/swirl)
