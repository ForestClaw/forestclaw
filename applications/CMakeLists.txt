#cmake_minimum_required(VERSION 3.13...3.20)

#project(ForestClawApps LANGUAGES C CXX Fortran)

#include(CTest)

#find_package(MPI COMPONENTS C Fortran REQUIRED)
#find_package(ZLIB REQUIRED)

# --- get installed ForestClaw package
#set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

#find_package(SC REQUIRED HINTS ${FORESTCLAW_ROOT} ${P4EST_ROOT} PATHS ${PROJECT_SOURCE_DIR}/../build)
#
#find_package(P4EST REQUIRED HINTS ${FORESTCLAW_ROOT} PATHS ${PROJECT_SOURCE_DIR}/../build)

#find_package(FORESTCLAW REQUIRED PATHS ${PROJECT_SOURCE_DIR}/../build)

# --- apps

add_subdirectory(clawpack/advection/2d/annulus)
add_subdirectory(clawpack/advection/2d/swirl)
