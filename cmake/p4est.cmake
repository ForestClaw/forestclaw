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
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${P4EST_ROOT} 
           -Dmpi:BOOL=${mpi} 
           -Dopenmp:BOOL=${openmp}
           -DSC_ROOT:PATH=${SC_ROOT}
           -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
BUILD_BYPRODUCTS ${P4EST_LIBRARIES}
DEPENDS SC-install
)
ExternalProject_Add_StepTargets(P4EST install)

# --- imported target

file(MAKE_DIRECTORY ${P4EST_INCLUDE_DIRS})
# avoid race condition

# this GLOBAL is required to be visible via other
# project's FetchContent of this project
add_library(P4EST::P4EST STATIC IMPORTED GLOBAL)
set_target_properties(P4EST::P4EST PROPERTIES 
  IMPORTED_LOCATION ${P4EST_LIBRARIES}
  INTERFACE_INCLUDE_DIRECTORIES ${P4EST_INCLUDE_DIRS}
  INTERFACE_LINK_LIBRARIES $<LINK_ONLY:SC::SC>
)

add_dependencies(P4EST::P4EST P4EST)
