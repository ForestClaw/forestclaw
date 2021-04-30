# provides imported targets P4EST::P4EST, ...
include(ExternalProject)
include(${CMAKE_CURRENT_LIST_DIR}/GitSubmodule.cmake)

set(p4est_external true CACHE BOOL "build p4est library" FORCE)

git_submodule("${PROJECT_SOURCE_DIR}/p4est")

# --- p4est externalProject
# this keeps project scopes totally separate, which avoids
# tricky to diagnose behaviors

if(NOT DEFINED P4EST_ROOT)
  set(P4EST_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

if(BUILD_SHARED_LIBS)
  set(P4EST_LIBRARIES ${P4EST_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}p4est${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(P4EST_LIBRARIES ${P4EST_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}p4est${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(P4EST_INCLUDE_DIRS ${P4EST_ROOT}/include)

ExternalProject_Add(P4EST
SOURCE_DIR ${PROJECT_SOURCE_DIR}/p4est
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${P4EST_ROOT} -Dmpi:BOOL=${mpi} -Dopenmp:BOOL=${openmp}
BUILD_BYPRODUCTS ${P4EST_LIBRARIES}
)

# --- imported target

file(MAKE_DIRECTORY ${P4EST_INCLUDE_DIRS})
# avoid race condition

# this GLOBAL is required to be visible via other
# project's FetchContent of this project
add_library(P4EST::P4EST INTERFACE IMPORTED GLOBAL)
target_include_directories(P4EST::P4EST INTERFACE "${P4EST_INCLUDE_DIRS}")
target_link_libraries(P4EST::P4EST INTERFACE "${P4EST_LIBRARIES}")

add_dependencies(P4EST::P4EST P4EST)
