# based on github.com/gemini3d/external
include_guard()

include(GNUInstallDirs)

include(${CMAKE_CURRENT_LIST_DIR}/ExtProj.cmake)

if(TARGET ZLIB::ZLIB)
  # if you don't want a find_package(ZLIB), must have project option to disable all find_package(ZLIB) calls
  # CMake cannot delete imported targets.
  return()
endif()

set(zlib_cmake_args
-DZLIB_COMPAT:BOOL=on
-DZLIB_ENABLE_TESTS:BOOL=off
-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=on
)
# CMAKE_POSITION_INDEPENDENT_CODE=on is needed for Zlib even when using static libs.

# before extproj for BUILD_BYPRODUCTS
if(BUILD_SHARED_LIBS)
  if(WIN32)
    set(ZLIB_LIBRARIES ${CMAKE_INSTALL_FULL_BINDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}zlib1${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(ZLIB_LIBRARIES ${CMAKE_INSTALL_FULL_LIBDIR}/${CMAKE_SHARED_LIBRARY_PREFIX}z${CMAKE_SHARED_LIBRARY_SUFFIX})
  endif()
else()
  if(MSVC)
    set(ZLIB_LIBRARIES ${CMAKE_INSTALL_FULL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}zlibstatic${CMAKE_STATIC_LIBRARY_SUFFIX})
  else()
    set(ZLIB_LIBRARIES ${CMAKE_INSTALL_FULL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}z${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()
endif()

extproj(zlib "${zlib_cmake_args}" "${ZLIB_LIBRARIES}" "")

# ZLIB::ZLIB

set(ZLIB_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include)

file(MAKE_DIRECTORY ${ZLIB_INCLUDE_DIRS})
# avoid race condition

add_library(ZLIB::ZLIB INTERFACE IMPORTED)
add_dependencies(ZLIB::ZLIB zlib)
target_link_libraries(ZLIB::ZLIB INTERFACE ${ZLIB_LIBRARIES})
target_include_directories(ZLIB::ZLIB INTERFACE ${ZLIB_INCLUDE_DIRS})
