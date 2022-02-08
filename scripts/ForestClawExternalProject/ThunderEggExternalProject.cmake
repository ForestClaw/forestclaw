#
# ThunderEggExternal_GetTarget
# -------------------------
#
# Define the imported CMake target for ThunderEgg
#
#     ThunderEggExternalProject_GetTarget(target
#                          [COMPONENTS arg1...]
#                          [OPTIONAL_COMPONENTS arg1...]
#     )
#
#
#   The options are:
#
#   COMPONENTS
#     list of components to build ThunderEgg with.
#     Current components are P4EST, FFTW, PETSC, LAPACK
#
#   OPTIONAL_COMPONENTS
#     list of optional components to build ThunderEgg with
#

function(ThunderEggExternalProject_GetTarget)
  set(one_value_args "")
  set(multi_value_args COMPONENTS OPTIONAL_COMPONENTS)
  cmake_parse_arguments(THUNDEREGG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN} )

  if(NOT DEFINED THUNDEREGG_ROOT)
    set(THUNDEREGG_ROOT ${CMAKE_INSTALL_PREFIX})
  endif()

  if(BUILD_SHARED_LIBS)
    set(THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ThunderEgg${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ThunderEgg${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()

  set(THUNDEREGG_INCLUDE_DIRS ${THUNDEREGG_ROOT}/include)

  find_package(MPI COMPONENTS CXX REQUIRED)

  set(THUNDEREGG_DEP_LIBS MPI::MPI_CXX)

  # --- optional dependent libraries
  if("PETSC" IN_LIST THUNDEREGG_COMPONENTS)
    find_package(PETSc REQUIRED)
  elseif("PETSC" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    find_package(PETSc)
  endif()

  if(TARGET PETSC::PETSC)
    list(APPEND THUNDEREGG_DEP_LIBS PETSC::PETSC)
  endif()


  if("FFTW" IN_LIST THUNDEREGG_COMPONENTS)
    find_package(FFTW REQUIRED)
  elseif("FFTW" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    find_package(FFTW)
  endif()

  if(TARGET FFTW::FFTW)
    list(APPEND THUNDEREGG_DEP_LIBS FFTW::FFTW)
  endif()

  if("P4EST" IN_LIST THUNDEREGG_COMPONENTS OR "P4EST" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    find_package(P4EST)
    find_package(SC)
  endif()

  if(TARGET P4EST::P4EST)
    list(APPEND THUNDEREGG_DEP_LIBS P4EST::P4EST)
  elseif("P4EST" IN_LIST THUNDEREGG_COMPONENTS OR "P4EST" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    if(BUILD_SHARED_LIBS)
      list(APPEND THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}p4est${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
      list(APPEND THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}p4est${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()
  endif()

  if(TARGET SC::SC)
    list(APPEND THUNDEREGG_DEP_LIBS SC::SC)
  elseif("P4EST" IN_LIST THUNDEREGG_COMPONENTS OR "P4EST" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    if(BUILD_SHARED_LIBS)
      list(APPEND THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}sc${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
      list(APPEND THUNDEREGG_LIBRARIES ${THUNDEREGG_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sc${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()
  endif()

  if("LAPACK" IN_LIST THUNDEREGG_COMPONENTS)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
  elseif("LAPACK" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    find_package(BLAS)
    find_package(LAPACK)
  endif()
  if(TARGET BLAS::BLAS)
    list(APPEND THUNDEREGG_DEP_LIBS BLAS::BLAS)
  endif()
  if(TARGET LAPACK::LAPACK)
    list(APPEND THUNDEREGG_DEP_LIBS LAPACK::LAPACK)
  endif()


  # --- imported target

  file(MAKE_DIRECTORY ${THUNDEREGG_INCLUDE_DIRS})
  # avoid race condition

  add_library(ThunderEgg::ThunderEgg STATIC IMPORTED GLOBAL)
  set_target_properties(ThunderEgg::ThunderEgg PROPERTIES 
    IMPORTED_LOCATION ${THUNDEREGG_LIBRARIES}
    INTERFACE_INCLUDE_DIRECTORIES ${THUNDEREGG_INCLUDE_DIRS}
    INTERFACE_LINK_LIBRARIES "${THUNDEREGG_DEP_LIBS}"
  )

endfunction(ThunderEggExternalProject_GetTarget)

#
# ThunderEggExternalProject
# -------------------------
#
# This module defines a function for using ThunderEgg as an external project
#
#     ThunderEggExternalProject(target
#                          TAG tag
#                          [REPOSITORY repo_url]
#                          [COMPONENTS arg1...]
#                          [OPTIONAL_COMPONENTS arg1...]
#     )
#
#
#   The options are:
#
#   TAG
#     Specifies the git tag from the repository to use
#
#   REPOSITORY
#     Specifies the url of git repository to use
#
#   COMPONENTS
#     list of components to build ThunderEgg with.
#     Current components are P4EST, FFTW, PETSC, LAPACK
#
#   OPTIONAL_COMPONENTS
#     list of optional components to build ThunderEgg with
#
function(ThunderEggExternalProject)
  include(ExternalProject)
  
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
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -Dpetsc=on -Dpetsc_required=on)
  elseif("PETSC" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -Dpetsc=on)
  else()
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -Dpetsc=off)
  endif()

  if(DEFINED PETSC_DIR)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DPETSC_DIR=${PETSC_DIR})
  endif()
  if(DEFINED PETSC_ARCH)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DPETSC_ARCH=${PETSC_ARCH})
  endif()


  if("FFTW" IN_LIST THUNDEREGG_COMPONENTS)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -Dfftw=on -Dfftw_required=on)
  elseif("FFTW" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -Dfftw=on)
  else()
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -Dfftw=off)
  endif()

  if(DEFINED FFTW_ROOT)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DFFTW_ROOT=${FFTW_ROOT})
  endif()

  if("P4EST" IN_LIST THUNDEREGG_COMPONENTS OR "P4EST" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -p4est=on)
  else()
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -p4est=off)
  endif()

  if(DEFINED P4EST_ROOT)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DP4EST_ROOT=${P4EST_ROOT})
  endif()

  if(DEFINED SC_ROOT)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -DSC_ROOT=${SC_ROOT})
  endif()

  if("LAPACK" IN_LIST THUNDEREGG_COMPONENTS)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -Dlapack=on -Dlapack_required=on)
  elseif("LAPACK" IN_LIST THUNDEREGG_OPTIONAL_COMPONENTS)
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -Dlapack=on)
  else()
    list(APPEND THUNDEREGG_CONFIGURE_ARGS -Dlapack=off)
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

  ThunderEggExternalProject_GetTarget(COMPONENTS ${THUNDEREGG_COMPONENTS} OPTIONAL_COMPONENTS ${THUNDEREGG_OPTIONAL_COMPONENTS})

  add_dependencies(ThunderEgg::ThunderEgg ThunderEgg)

endfunction(ThunderEggExternalProject)