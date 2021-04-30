cmake_minimum_required(VERSION 3.13...3.20)

project(ForestClawApps LANGUAGES C CXX Fortran)

include(CTest)

find_package(ZLIB REQUIRED)
find_package(MPI COMPONENTS C CXX Fortran REQUIRED)


# --- get installed ForestClaw package
find_package(FORESTCLAW REQUIRED)

# build convenience imported target to assure proper linking
add_library(FORESTCLAW::FORESTCLAW INTERFACE IMPORTED)
target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE
FORESTCLAW::clawpack4.6 FORESTCLAW::clawpack5 FORESTCLAW::clawpatch
P4EST::P4EST SC::SC MPI::MPI_Fortran MPI::MPI_CXX MPI::MPI_C ZLIB::ZLIB)

# --- apps

add_subdirectory(clawpack/advection/2d/swirl)
