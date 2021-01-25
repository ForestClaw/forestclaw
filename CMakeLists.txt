cmake_minimum_required(VERSION 3.13)

include(cmake/policy.cmake)

project(FORESTCLAW
  LANGUAGES C Fortran
  VERSION 0.1.4880)

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

if(NOT p4est_external)
  find_package(P4EST)
endif()
if(NOT TARGET P4EST::P4EST)
  include(cmake/p4est.cmake)
endif()

# --- ForestClaw

include(cmake/config.cmake)
include(cmake/compilers.cmake)

add_library(forestclaw)
# source and links defined in src/CMakeLists.txt

# imported target, for use from FetchContent
add_library(FORESTCLAW::FORESTCLAW INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE forestclaw)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ ${PROJECT_BINARY_DIR}/include/
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")

install(TARGETS forestclaw
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


add_subdirectory(src)

include(cmake/install.cmake)
