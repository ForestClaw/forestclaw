include(CheckIncludeFile)
include(CheckIncludeFiles)
include(CheckSymbolExists)

# --- generate fclaw_config.h

set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)

if(MPI_FOUND)
  set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_C)
  set(FCLAW_CC \"${MPI_C_COMPILER}\")
  set(FCLAW_CPP ${MPI_C_COMPILER})
  set(FCLAW_CXX \"${MPI_CXX_COMPILER}\")
  SET(FCLAW_F77 \"${MPI_Fortran_COMPILER}\")
else()
  set(FCLAW_CC \"${CMAKE_C_COMPILER}\")
  set(FCLAW_CPP ${CMAKE_C_COMPILER})
  set(FCLAW_CXX \"${CMAKE_CXX_COMPILER}\")
  SET(FCLAW_F77 \"${CMAKE_Fortran_COMPILER}\")
endif()

string(APPEND FCLAW_CPP " -E")
set(FCLAW_CPP \"${FCLAW_CPP}\")

set(FCLAW_CFLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_OPTIONS}")
set(FCLAW_CFLAGS \"${FCLAW_CFLAGS}\")

set(FCLAW_CPPFLAGS \"\")

set(FCLAW_FFLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_OPTIONS}")
set(FCLAW_FFLAGS \"${FCLAW_FFLAGS}\")

set(FCLAW_FLIBS \"${MPI_Fortran_LIBRARIES}\")

set(FCLAW_LDFLAGS \"${MPI_C_LINK_FLAGS}\")
set(FCLAW_LIBS \"${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${ZLIB_LIBRARIES} m\")

set(FCLAW_ENABLE_MEMALIGN 1)

if(MPI_FOUND)
  set(FCLAW_ENABLE_MPI ${MPI_FOUND})
  set(FCLAW_ENABLE_MPIIO 1)
endif(MPI_FOUND)

# check_symbol_exists(sqrt math.h FCLAW_NONEED_M)
# if(NOT FCLAW_NONEED_M)
#   set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} m)
#   check_symbol_exists(sqrt math.h FCLAW_NEED_M)
# endif()

check_include_file(fenv.h FCLAW_HAVE_FENV_H)
if(FCLAW_HAVE_FENV_H)
  set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
  if(NOT WIN32)
    set(CMAKE_REQUIRED_LIBRARIES m)
  endif()
  check_symbol_exists(feenableexcept fenv.h FCLAW_HAVE_FEENABLEEXCEPT)
  set(CMAKE_REQUIRED_LIBRARIES)
endif()

check_include_file(signal.h FCLAW_HAVE_SIGNAL_H)

check_include_file(unistd.h FCLAW_HAVE_UNISTD_H)

set(FCLAW_PACKAGE \"${PROJECT_NAME}\")

if(TARGET HDF5::HDF5)
  set(FCLAW_ENABLE_HDF5 1)
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/fclaw_config.h.in ${PROJECT_BINARY_DIR}/include/fclaw_config.h)
install(FILES ${PROJECT_BINARY_DIR}/include/fclaw_config.h TYPE INCLUDE)
configure_file(${CMAKE_CURRENT_LIST_DIR}/test_config.h.in ${PROJECT_BINARY_DIR}/test/test_config.h)

set(top_builddir ${PROJECT_BINARY_DIR})
set(top_srcdir ${PROJECT_SOURCE_DIR})
configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile)
