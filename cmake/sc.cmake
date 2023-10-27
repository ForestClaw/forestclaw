# provides imported target SC::SC
include(ExternalProject)
include(GNUInstallDirs)

include(${CMAKE_CURRENT_LIST_DIR}/GitSubmodule.cmake)

git_submodule("${PROJECT_SOURCE_DIR}/sc")

# --- libsc externalProject
# this keeps project scopes totally separate, which avoids
# tricky to diagnose behaviors

if(BUILD_SHARED_LIBS)
  if(WIN32)
    set(SC_LIBRARIES ${CMAKE_INSTALL_PREFIX}/bin/${CMAKE_SHARED_LIBRARY_PREFIX}sc${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(SC_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}sc${CMAKE_SHARED_LIBRARY_SUFFIX})
  endif()
else()
  set(SC_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sc${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(SC_INCLUDE_DIRS ${CMAKE_INSTALL_FULL_INCLUDEDIR})

set(cmake_sc_args
-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
-Dmpi:BOOL=${mpi}
-Dopenmp:BOOL=${openmp}
-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
)

ExternalProject_Add(SC
SOURCE_DIR ${PROJECT_SOURCE_DIR}/sc
CMAKE_ARGS ${cmake_sc_args}
BUILD_BYPRODUCTS ${SC_LIBRARIES}
TEST_COMMAND ""
CONFIGURE_HANDLED_BY_BUILD ON
INACTIVITY_TIMEOUT 60
USES_TERMINAL_DOWNLOAD true
USES_TERMINAL_UPDATE true
USES_TERMINAL_PATCH true
USES_TERMINAL_CONFIGURE true
USES_TERMINAL_BUILD true
USES_TERMINAL_INSTALL true
USES_TERMINAL_TEST true
)

# --- imported target

file(MAKE_DIRECTORY ${SC_INCLUDE_DIRS})
# avoid race condition

# this GLOBAL is required to be visible via other
# project's FetchContent of this project
add_library(SC::SC INTERFACE IMPORTED GLOBAL)
target_link_libraries(SC::SC INTERFACE ${SC_LIBRARIES} ZLIB::ZLIB)
target_include_directories(SC::SC INTERFACE ${SC_INCLUDE_DIRS})

add_dependencies(SC::SC SC)
