# provides imported target SC::SC
include(ExternalProject)
include(${CMAKE_CURRENT_LIST_DIR}/GitSubmodule.cmake)

set(sc_external true CACHE BOOL "build sc library" FORCE)

git_submodule("${PROJECT_SOURCE_DIR}/sc")

# --- libsc externalProject
# this keeps project scopes totally separate, which avoids
# tricky to diagnose behaviors

if(NOT DEFINED SC_ROOT)
  set(SC_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

if(BUILD_SHARED_LIBS)
  if(WIN32)
    set(SC_LIBRARIES ${SC_ROOT}/bin/${CMAKE_SHARED_LIBRARY_PREFIX}sc${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(SC_LIBRARIES ${SC_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}sc${CMAKE_SHARED_LIBRARY_SUFFIX})
  endif()
else()
  set(SC_LIBRARIES ${SC_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sc${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(SC_INCLUDE_DIRS ${SC_ROOT}/include)

set(cmake_sc_args
-DCMAKE_INSTALL_PREFIX:PATH=${SC_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
-Dmpi:BOOL=${mpi}
-Dopenmp:BOOL=${openmp}
)

ExternalProject_Add(SC
SOURCE_DIR ${PROJECT_SOURCE_DIR}/sc
CMAKE_ARGS ${cmake_sc_args}
BUILD_BYPRODUCTS ${SC_LIBRARIES}
)
ExternalProject_Add_StepTargets(SC install)

# --- imported target

file(MAKE_DIRECTORY ${SC_INCLUDE_DIRS})
# avoid race condition

# this GLOBAL is required to be visible via other
# project's FetchContent of this project
add_library(SC::SC STATIC IMPORTED GLOBAL)
set_target_properties(SC::SC PROPERTIES 
  IMPORTED_LOCATION ${SC_LIBRARIES}
  INTERFACE_INCLUDE_DIRECTORIES ${SC_INCLUDE_DIRS}
  INTERFACE_LINK_LIBRARIES $<LINK_ONLY:ZLIB::ZLIB>
)

add_dependencies(SC::SC SC)
