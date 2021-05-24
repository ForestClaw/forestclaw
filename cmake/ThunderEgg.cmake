# provides imported targets THUNDEREGG::THUNDEREGG, ...
include(ExternalProject)
include(${CMAKE_CURRENT_LIST_DIR}/GitSubmodule.cmake)

set(thunderegg_external true CACHE BOOL "build thunderegg library" FORCE)

git_submodule("${PROJECT_SOURCE_DIR}/ThunderEgg")

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

ExternalProject_Add(THUNDEREGG
SOURCE_DIR ${PROJECT_SOURCE_DIR}/thunderegg
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${THUNDEREGG_ROOT} -Dmpi:BOOL=${mpi} -Dopenmp:BOOL=${openmp}
BUILD_BYPRODUCTS ${THUNDEREGG_LIBRARIES}
)

# --- imported target

file(MAKE_DIRECTORY ${THUNDEREGG_INCLUDE_DIRS})
# avoid race condition

# this GLOBAL is required to be visible via other
# project's FetchContent of this project
add_library(THUNDEREGG::THUNDEREGG INTERFACE IMPORTED GLOBAL)
target_include_directories(THUNDEREGG::THUNDEREGG INTERFACE "${THUNDEREGG_INCLUDE_DIRS}")
target_link_libraries(THUNDEREGG::THUNDEREGG INTERFACE "${THUNDEREGG_LIBRARIES}")

add_dependencies(THUNDEREGG::THUNDEREGG THUNDEREGG)
