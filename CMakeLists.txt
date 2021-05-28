cmake_minimum_required(VERSION 3.13...3.20)

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
  find_package(MPI COMPONENTS C REQUIRED)
endif()

find_package(LAPACK)
find_package(BLAS)
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

if(NOT thunderegg_external)
  find_package(ThunderEgg HINTS ${THUNDEREGG_ROOT})
endif()
if(NOT TARGET THUNDEREGG::THUNDEREGG)
  include(cmake/ThunderEgg.cmake)
endif()

# --- configure

include(cmake/config.cmake)
include(cmake/compilers.cmake)

# --- test

add_subdirectory(test)

# --- ForestClaw

add_library(forestclaw)
set_target_properties(forestclaw PROPERTIES EXPORT_NAME FORESTCLAW)
target_link_libraries(forestclaw PUBLIC SC::SC P4EST::P4EST)
target_link_libraries(forestclaw PRIVATE ZLIB::ZLIB)
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

if(clawpack)

  # --- Clawpatch

  add_library(clawpatch)
  target_link_libraries(clawpatch PUBLIC forestclaw)
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

  # imported target, for use from FetchContent
  add_library(FORESTCLAW::CLAWPATCH INTERFACE IMPORTED GLOBAL)
  target_link_libraries(FORESTCLAW::CLAWPATCH INTERFACE clawpatch)


  # --- Clawpack 4.6

  add_library(clawpack4.6)
  target_link_libraries(clawpack4.6 PUBLIC clawpatch)
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

  # imported target, for use from FetchContent
  add_library(FORESTCLAW::CLAWPACK4.6 INTERFACE IMPORTED GLOBAL)
  target_link_libraries(FORESTCLAW::CLAWPACK4.6 INTERFACE clawpack4.6)

  # --- Clawpack 5

  add_library(clawpack5)
  target_link_libraries(clawpack5 PUBLIC clawpatch)
  target_include_directories(clawpack5
    INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/solvers/fc2d_clawpack5>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include/solvers/fc2d_clawpack5>
    $<INSTALL_INTERFACE:include>)
  if(mpi)
    target_link_libraries(clawpack5 PRIVATE MPI::MPI_C)
  endif(mpi)


  install(DIRECTORY
    ${PROJECT_SOURCE_DIR}/src/solvers/fc2d_clawpack5
    DESTINATION include
    FILES_MATCHING PATTERN "*.h")

  install(TARGETS clawpack5
    EXPORT ${PROJECT_NAME}Targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)

  add_subdirectory(src/solvers/fc2d_clawpack5)

  # imported target, for use from FetchContent
  add_library(FORESTCLAW::CLAWPACK5 INTERFACE IMPORTED GLOBAL)
  target_link_libraries(FORESTCLAW::CLAWPACK5 INTERFACE clawpack5)

endif(clawpack)

# --- geoclaw

if(geoclaw)

  add_library(geoclaw)
  set_target_properties(geoclaw PROPERTIES EXPORT_NAME GEOCLAW)
  target_link_libraries(geoclaw PUBLIC forestclaw clawpatch)

  add_subdirectory(src/solvers/fc2d_geoclaw)

  # imported target, for use from FetchContent
  add_library(FORESTCLAW::GEOCLAW INTERFACE IMPORTED GLOBAL)
  target_link_libraries(FORESTCLAW::GEOCLAW INTERFACE clawpatch)

endif(geoclaw)

# --- ThunderEgg
if(thunderegg)

  add_library(thunderegg)
  set_target_properties(thunderegg PROPERTIES EXPORT_NAME THUNDEREGG)
  target_link_libraries(thunderegg PUBLIC forestclaw clawpatch ThunderEgg::ThunderEgg)
  if(mpi)
    target_link_libraries(thunderegg PRIVATE MPI::MPI_C)
  endif(mpi)

  add_subdirectory(src/solvers/fc2d_thunderegg)

  # imported target, for use from FetchContent
  add_library(FORESTCLAW::THUNDEREGG INTERFACE IMPORTED GLOBAL)
  target_link_libraries(FORESTCLAW::THUNDEREGG INTERFACE thunderegg)

endif(thunderegg)
# --- package

include(cmake/pkgconf.cmake)
include(cmake/install.cmake)

# --- examples
add_subdirectory(applications)
