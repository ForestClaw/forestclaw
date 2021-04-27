include(FetchContent)

# provides P4EST::P4EST

set(p4est_external true CACHE BOOL "build p4est library" FORCE)

include(${CMAKE_CURRENT_LIST_DIR}/sc.cmake)

# FIXME: switch to p4est/ submodule once our CMake PR is merged in p4est.
FetchContent_Declare(P4EST
  GIT_REPOSITORY https://github.com/scivision/p4est.git
  GIT_TAG feature-cmake
  CMAKE_ARGS -Dmpi=${mpi})

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  FetchContent_MakeAvailable(P4EST)
elseif(NOT p4est_POPULATED)
  FetchContent_Populate(P4EST)
  add_subdirectory(${p4est_SOURCE_DIR} ${p4est_BINARY_DIR})
endif()

target_link_libraries(P4EST::P4EST INTERFACE SC::SC)
