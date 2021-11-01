# provides imported target P4EST::P4EST

function(ThunderEggExternalProject)
  include(ExternalProject)
  
  set(options REQUIRED)
  set(one_value_args REPOSITORY TAG)
  set(multi_value_args COMPONENTS OPTIONAL_COMPONENTS)
  cmake_parse_arguments(THUNDEREGG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN} )

  if(NOT DEFINED THUNDEREGG_ROOT)
    set(THUNDEREGG_ROOT ${CMAKE_INSTALL_PREFIX})
  endif()

  if(NOT DEFINED THUNDEREGG_REPOSITORY)
    set(THUNDEREGG_REPOSITORY "https://github.com/thunderegg/thunderegg")
  endif()

  if(NOT DEFINED THUNDEREGG_TAG)
    message(FATAL_ERROR "TAG for ThunderEgg External Project needs to be set")
  endif()

  if(BUILD_SHARED_LIBS)
    set(THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ThunderEgg${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ThunderEgg${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()

  set(THUNDEREGG_INCLUDE_DIRS ${THUNDEREGG_ROOT}/include)

  find_package(MPI COMPONENTS CXX REQUIRED)


  set(THUNDEREGG_CONFIGURE_ARGS -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=${THUNDEREGG_ROOT})
  set(THUNDEREGG_DEP_LIBS MPI::MPI_CXX)

  # --- optional dependent libraries
  if("PETSC" IN_LIST THUNDEREGG_COMPONENTS)
    find_package(PETSc REQUIRED)
  else()
    find_package(PETSc)
  endif()
  if(DEFINED PETSC_DIR)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DPETSC_DIR=${PETSC_DIR})
  endif()
  if(DEFINED PETSC_ARCH)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DPETSC_ARCH=${PETSC_ARCH})
  endif()
  if(TARGET PETSC::PETSC)
    list(APPEND THUNDEREGG_DEP_LIBS PETSC::PETSC)
  endif()


  if("FFTW" IN_LIST THUNDEREGG_COMPONENTS)
    find_package(FFTW REQUIRED)
  else()
    find_package(FFTW)
  endif()
  if(DEFINED FFTW_ROOT)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DFFTW_ROOT=${FFTW_ROOT})
  endif()
  if(TARGET FFTW::FFTW)
    list(APPEND THUNDEREGG_DEP_LIBS FFTW::FFTW)
  endif()

  find_package(P4EST)
  if(DEFINED P4EST_ROOT)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DP4EST_ROOT=${P4EST_ROOT})
  endif()
  if(TARGET P4EST::P4EST)
    list(APPEND THUNDEREGG_DEP_LIBS P4EST::P4EST)
  else()
    if(BUILD_SHARED_LIBS)
      list(APPEND THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}p4est${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
      list(APPEND THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}p4est${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()
  endif()

  find_package(SC)
  if(DEFINED SC_ROOT)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DSC_ROOT=${SC_ROOT})
  endif()
  if(TARGET SC::SC)
    list(APPEND THUNDEREGG_DEP_LIBS SC::SC)
  else()
    if(BUILD_SHARED_LIBS)
      list(APPEND THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}sc${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
      list(APPEND THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sc${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()
  endif()

  if("LAPACK" IN_LIST THUNDEREGG_COMPONENTS)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
  else()
    find_package(BLAS)
    find_package(LAPACK)
  endif()
  if(TARGET BLAS::BLAS)
    list(APPEND THUNDEREGG_DEP_LIBS BLAS::BLAS)
  endif()
  if(TARGET LAPACK::LAPACK)
    list(APPEND THUNDEREGG_DEP_LIBS LAPACK::LAPACK)
  endif()


  ExternalProject_Add(ThunderEgg
    GIT_REPOSITORY   ${THUNDEREGG_REPOSITORY}
    GIT_TAG          ${THUNDEREGG_TAG}
    CMAKE_ARGS       ${THUNDEREGG_CONFIGURE_ARGS}
    BUILD_BYPRODUCTS ${THUNDEREGG_LIBRARIES}
    DEPENDS          ${THUNDEREGG_EXTERNAL_PROJECT_DEPENDS}
  )

  # --- imported target

  file(MAKE_DIRECTORY ${THUNDEREGG_INCLUDE_DIRS})
  # avoid race condition

  # this GLOBAL is required to be visible via other
  # project's FetchContent of this project
  add_library(ThunderEgg::ThunderEgg STATIC IMPORTED GLOBAL)
  set_target_properties(ThunderEgg::ThunderEgg PROPERTIES 
    IMPORTED_LOCATION ${THUNDEREGG_LIBRARIES}
    INTERFACE_INCLUDE_DIRECTORIES ${THUNDEREGG_INCLUDE_DIRS}
    INTERFACE_LINK_LIBRARIES "${THUNDEREGG_DEP_LIBS}"
  )

  add_dependencies(ThunderEgg::ThunderEgg ThunderEgg)

endfunction(ThunderEggExternalProject)