# provides imported targets THUNDEREGG::THUNDEREGG, ...
include(ExternalProject)

set(thunderegg_external true CACHE BOOL "build thunderegg library" FORCE)

# --- thunderegg externalProject
# this keeps project scopes totally separate, which avoids
# tricky to diagnose behaviors

if(NOT DEFINED THUNDEREGG_ROOT)
  set(THUNDEREGG_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

if(BUILD_SHARED_LIBS)
  set(THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ThunderEgg${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ThunderEgg${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(THUNDEREGG_INCLUDE_DIRS ${THUNDEREGG_ROOT}/include)

if(p4est_external)
  set(THUNDEREGG_EXTERNAL_PROJECT_DEPENDS "P4EST-install")
endif()

ExternalProject_Add(ThunderEgg
GIT_REPOSITORY https://github.com/thunderegg/thunderegg
GIT_TAG        v1.0.3
CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX:PATH=${THUNDEREGG_ROOT} -Dmpi:BOOL=${mpi} -Dopenmp:BOOL=${openmp} -Ddisable_petsc:BOOL=true -DP4EST_ROOT=${P4EST_ROOT} -DBUILD_TESTING:BOOL=off
BUILD_BYPRODUCTS ${THUNDEREGG_LIBRARIES}
DEPENDS ${THUNDEREGG_EXTERNAL_PROJECT_DEPENDS}
)

# --- required libraries
find_package(FFTW REQUIRED)
find_package(BLAS)
find_package(LAPACK)
find_package(MPI COMPONENTS CXX)

# --- imported target

file(MAKE_DIRECTORY ${THUNDEREGG_INCLUDE_DIRS})
# avoid race condition

# this GLOBAL is required to be visible via other
# project's FetchContent of this project
add_library(ThunderEgg::ThunderEgg STATIC IMPORTED GLOBAL)
target_include_directories(ThunderEgg::ThunderEgg INTERFACE "${THUNDEREGG_INCLUDE_DIRS}")
target_link_libraries(ThunderEgg::ThunderEgg INTERFACE "${THUNDEREGG_LIBRARIES}" P4EST::P4EST)
set_target_properties(ThunderEgg::ThunderEgg PROPERTIES 
  IMPORTED_LOCATION ${THUNDEREGG_LIBRARIES}
  INTERFACE_INCLUDE_DIRECTORIES ${THUNDEREGG_INCLUDE_DIRS}
  INTERFACE_LINK_LIBRARIES "FFTW::FFTW;P4EST::P4EST;SC::SC;MPI::MPI_CXX"
)
if(TARGET BLAS::BLAS AND TARGET LAPACK::LAPACK)
  target_link_libraries(ThunderEgg::ThunderEgg INTERFACE BLAS::BLAS LAPACK::LAPACK)
endif()

add_dependencies(ThunderEgg::ThunderEgg ThunderEgg)
