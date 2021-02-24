cmake_minimum_required(VERSION 3.13...3.20)

project(FORESTCLAW
  LANGUAGES C CXX Fortran
  VERSION 0.1.4880)

message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}")

include(CTest)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# --- user options

include(cmake/options.cmake)

# --- external libs

if(mpi)
  find_package(MPI COMPONENTS C REQUIRED)
endif()

find_package(LAPACK)
find_package(BLAS)
find_package(ZLIB)

# --- p4est, sc

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

if(NOT sc_external)
  find_package(SC HINTS ${P4EST_ROOT})
endif()
if(NOT TARGET SC::SC)
  include(cmake/sc.cmake)
endif()

if(NOT p4est_external)
  find_package(P4EST)
endif()
if(NOT TARGET P4EST::P4EST)
  include(cmake/p4est.cmake)
endif()

# --- configure

include(cmake/config.cmake)
include(cmake/compilers.cmake)

# --- ForestClaw

add_library(forestclaw)
set_target_properties(forestclaw PROPERTIES EXPORT_NAME FORESTCLAW)
target_include_directories(forestclaw
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)

# imported target, for use from FetchContent
add_library(FORESTCLAW::FORESTCLAW INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE forestclaw)

install(DIRECTORY
  ${PROJECT_SOURCE_DIR}/src/
  ${PROJECT_BINARY_DIR}/include/
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")

install(TARGETS forestclaw
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

add_subdirectory(src)

# --- Clawpatch

add_library(clawpatch)
set_target_properties(clawpatch PROPERTIES EXPORT_NAME CLAWPATCH)
target_link_libraries(clawpatch PRIVATE forestclaw P4EST::P4EST SC::SC)
target_include_directories(clawpatch
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/patches/clawpatch>
  $<INSTALL_INTERFACE:include/clawpatch>
  $<INSTALL_INTERFACE:include>)

install(DIRECTORY
  ${PROJECT_SOURCE_DIR}/src/patches/clawpatch
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")

install(TARGETS clawpatch
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

add_subdirectory(src/patches/clawpatch)

# --- Clawpack 4.6

add_library(clawpack4.6)
set_target_properties(clawpack4.6 PROPERTIES EXPORT_NAME CLAWPACK4.6)
target_link_libraries(clawpack4.6 PRIVATE clawpatch P4EST::P4EST SC::SC)
target_include_directories(clawpack4.6
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/solvers/fc2d_clawpack4.6>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include/solvers/fc2d_clawpack4.6>
  $<INSTALL_INTERFACE:include>)

install(DIRECTORY
  ${PROJECT_SOURCE_DIR}/src/solvers/fc2d_clawpack4.6
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")

install(TARGETS clawpack4.6
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

add_subdirectory(src/solvers/fc2d_clawpack4.6)

# --- Clawpack 5

add_library(clawpack5)
set_target_properties(clawpack5 PROPERTIES EXPORT_NAME CLAWPACK5)
target_link_libraries(clawpack5 PRIVATE clawpatch P4EST::P4EST SC::SC)
target_include_directories(clawpack5
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/solvers/fc2d_clawpack5>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include/solvers/fc2d_clawpack5>
  $<INSTALL_INTERFACE:include>)

install(DIRECTORY
  ${PROJECT_SOURCE_DIR}/src/solvers/fc2d_clawpack5
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")

install(TARGETS clawpack5
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

add_subdirectory(src/solvers/fc2d_clawpack5)

# --- geoclaw

if(geoclaw)

add_library(geoclaw)
target_link_libraries(geoclaw PRIVATE forestclaw clawpatch P4EST::P4EST SC::SC)

add_subdirectory(src/solvers/fc2d_geoclaw)


endif(geoclaw)

# --- package

include(cmake/install.cmake)
