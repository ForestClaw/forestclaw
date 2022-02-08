# FindFFTW.cmake
# ---------------
#
# Find fftw library
#
# Result Variables
# ----------------
#
# This module defines the following variables::
#
#   FFTW_FOUND
#   FFTW_INCLUDE_DIRS   - include directories for p4est
#   FFTW_LIBRARIES      - link against this library to use p4est
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#   FFTW::FFTW

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if(PKG_CONFIG_FOUND)
  pkg_check_modules( PKG_FFTW QUIET "fftw3" )
endif(PKG_CONFIG_FOUND)

if(NOT FFTW_FOUND)

  find_path (FFTW_INCLUDE_DIR
    NAMES fftw3.h
    DOC "fftw3 header")
  
  find_library (FFTW_LIBRARY
    NAMES fftw3
    DOC "fftw3 library")
  
  if(FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
    set(FFTW_FFTW_FOUND true)
  endif()
  
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args (FFTW
    REQUIRED_VARS FFTW_LIBRARY FFTW_INCLUDE_DIR
    HANDLE_COMPONENTS)
  
endif(NOT FFTW_FOUND)

if(FFTW_FOUND)

  set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
  set(FFTW_LIBRARIES ${FFTW_LIBRARY})
  
  if(NOT TARGET FFTW::FFTW)
      add_library(FFTW::FFTW INTERFACE IMPORTED)
      set_target_properties(FFTW::FFTW PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
          INTERFACE_LINK_LIBRARIES "${FFTW_LIBRARY}"
      )
  endif(NOT TARGET FFTW::FFTW)

endif(FFTW_FOUND)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY)