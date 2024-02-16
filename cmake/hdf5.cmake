# based on github.com/gemini3d/external
include(GNUInstallDirs)

include(${CMAKE_CURRENT_LIST_DIR}/ExtProj.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/zlib.cmake)

# pass MPI hints to HDF5, which needs it as a CMake variable--HDF5 project doesn't detect ENV{MPI_ROOT}
if(NOT MPI_ROOT AND DEFINED ENV{MPI_ROOT})
  set(MPI_ROOT $ENV{MPI_ROOT})
endif()
#TODO remove this
set(hdf5_parallel true)

set(hdf5_cmake_args
-DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
-DZLIB_USE_EXTERNAL:BOOL=OFF
-DZLIB_ROOT:PATH=${CMAKE_INSTALL_PREFIX}
-DHDF5_GENERATE_HEADERS:BOOL=false
-DHDF5_DISABLE_COMPILER_WARNINGS:BOOL=true
-DBUILD_STATIC_LIBS:BOOL=$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>
-DHDF5_BUILD_FORTRAN:BOOL=true
-DHDF5_BUILD_CPP_LIB:BOOL=false
-DHDF5_BUILD_EXAMPLES:BOOL=false
-DHDF5_BUILD_TOOLS:BOOL=$<NOT:$<BOOL:${hdf5_parallel}>>
-DHDF5_ENABLE_PARALLEL:BOOL=$<BOOL:${hdf5_parallel}>
)
# https://github.com/HDFGroup/hdf5/issues/818  for broken ph5diff in HDF5_BUILD_TOOLS
if(MPI_ROOT)
  list(APPEND hdf5_cmake_args -DMPI_ROOT:PATH=${MPI_ROOT})
endif()

# NOTE: don't use CMAKE_INSTALL_FULL_LIBDIR as HDF5 1.10 doesn't use that lib/lib64 convention
set(HDF5_LIBRARIES)
foreach(_name IN ITEMS hdf5_hl_fortran hdf5_hl_f90cstub hdf5_fortran hdf5_f90cstub hdf5_hl hdf5)
  if(BUILD_SHARED_LIBS)
    if(WIN32)
      list(APPEND HDF5_LIBRARIES ${CMAKE_INSTALL_FULL_BINDIR}/lib${_name}${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
      list(APPEND HDF5_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/lib${_name}${CMAKE_SHARED_LIBRARY_SUFFIX})
    endif()
  else()
    list(APPEND HDF5_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/lib${_name}${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()
endforeach()

set(hdf5_req)
if(NOT TARGET ZLIB::ZLIB)
  set(hdf5_req zlib)
endif()

extproj(hdf5 "${hdf5_cmake_args}" "${HDF5_LIBRARIES}" "${hdf5_req}")

# HDF5::HDF5

set(HDF5_INCLUDE_DIRS ${CMAKE_INSTALL_FULL_INCLUDEDIR})

file(MAKE_DIRECTORY ${HDF5_INCLUDE_DIRS})
# avoid race condition

add_library(HDF5::HDF5 INTERFACE IMPORTED)
add_dependencies(HDF5::HDF5 hdf5)
target_include_directories(HDF5::HDF5 INTERFACE "${HDF5_INCLUDE_DIRS}")
target_link_libraries(HDF5::HDF5 INTERFACE "${HDF5_LIBRARIES}")

# --- HDF5 parallel compression support
# this could be improved by making it an ExternalProject post-build step instead of assumptions made here
if(hdf5_parallel)
  if(MPI_C_VERSION VERSION_GREATER_EQUAL 3)
    message(STATUS "HDF5-MPI: MPI-3 available, assuming HDF5 parallel compression enabled")
    set(hdf5_parallel_compression true)
  else()
    message(STATUS "HDF5-MPI: MPI-3 NOT available => HDF5 parallel compression disabled")
    set(hdf5_parallel_compression false)
  endif()
endif()

# --- external deps


find_package(Threads)

target_link_libraries(HDF5::HDF5 INTERFACE
ZLIB::ZLIB
${CMAKE_THREAD_LIBS_INIT}
${CMAKE_DL_LIBS}
$<$<BOOL:${UNIX}>:m>
)
