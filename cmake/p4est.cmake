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

set(FCLAW_CONFIGURE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${P4EST_ROOT} -Dmpi:BOOL=${mpi} -Dopenmp:BOOL=${openmp} -DSC_ROOT=${SC_ROOT})

# --- compiler flags
list(APPEND FCLAW_CONFIGURE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS})
list(APPEND FCLAW_CONFIGURE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS})
list(APPEND FCLAW_CONFIGURE_ARGS -DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS})

# --- propagate flags to subprojects
GET_CMAKE_PROPERTY(CACHE_VARS CACHE_VARIABLES)
FOREACH(CACHE_VAR ${CACHE_VARS})
  GET_PROPERTY(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
  IF(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line.")
    GET_PROPERTY(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    IF(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
      SET(CACHE_VAR_TYPE)
    ELSE(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
      SET(CACHE_VAR_TYPE :${CACHE_VAR_TYPE})
    ENDIF()
    list(APPEND FCLAW_CONFIGURE_ARGS -D${CACHE_VAR}${CACHE_VAR_TYPE}=\"${${CACHE_VAR}}\")
  ENDIF()
ENDFOREACH(CACHE_VAR ${CACHE_VARS})

ExternalProject_Add(P4EST
SOURCE_DIR ${PROJECT_SOURCE_DIR}/p4est
CMAKE_ARGS ${P4EST_CONFIGURE_ARGS}
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
