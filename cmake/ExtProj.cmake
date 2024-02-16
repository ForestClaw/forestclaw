# based on github.com/gemini3d/external

include_guard()

include(ExternalProject)

if(CMAKE_VERSION VERSION_LESS 3.19)
  message(FATAL_ERROR "CMake >= 3.19 required to use JSON")
endif()

file(READ ${CMAKE_CURRENT_LIST_DIR}/libraries.json json)


function(extproj name cmake_args byproducts depends)

# PREPEND so that user arguments can override these defaults
list(PREPEND cmake_args
-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
-DCMAKE_PREFIX_PATH:PATH=${CMAKE_INSTALL_PREFIX}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
-DBUILD_TESTING:BOOL=false
-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
)
if(CMAKE_TOOLCHAIN_FILE)
  list(APPEND cmake_args -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE})
endif()

# builds each project in parallel, without needing to build all projects simultaneously in parallel.
# this greatly aids debugging while maintaining speed of build overall.
if(CMAKE_GENERATOR MATCHES "Makefiles" AND NOT DEFINED ENV{CMAKE_BUILD_PARALLEL_LEVEL})
  cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
endif()
set(build_parallel ${CMAKE_COMMAND} --build <BINARY_DIR> --parallel ${Ncpu})

string(JSON url GET ${json} ${name} url)
if(url MATCHES ".git$")

  string(JSON tag GET ${json} ${name} tag)

  set(download_args
  GIT_REPOSITORY ${url}
  GIT_TAG ${tag}
  GIT_SHALLOW true
  GIT_PROGRESS true
  )

else()

  set(download_args
  URL ${url}
  )

endif()


# BUILD_BYPRODUCTS is vital for Ninja, else build-time error "missing and no known rule to make it"

ExternalProject_Add(${name}
${download_args}
BUILD_COMMAND ${build_parallel}
TEST_COMMAND ""
CMAKE_ARGS ${cmake_args}
DEPENDS ${depends}
BUILD_BYPRODUCTS ${byproducts}
INACTIVITY_TIMEOUT 60
CONFIGURE_HANDLED_BY_BUILD true
USES_TERMINAL_DOWNLOAD true
USES_TERMINAL_UPDATE true
USES_TERMINAL_BUILD true
USES_TERMINAL_INSTALL true
USES_TERMINAL_TEST true
)

endfunction(extproj)
