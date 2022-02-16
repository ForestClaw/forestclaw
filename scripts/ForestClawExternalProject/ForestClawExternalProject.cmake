#
# ForestClawExternal_GetTarget
# ----------------------------
#
# Define the imported CMake target for ForestClaw
#
#     ForestClawExternalProject_GetTarget(target
#                                        [COMPONENTS arg1...]
#                                        )
#
#
#   The options are:
#
#   COMPONENTS
#     list of components to build ForestClaw with.
#     Current components are P4EST, FFTW, PETSC, LAPACK
#
function(ForestClawExternalProject_GetTarget)
  set(one_value_args "")
  set(multi_value_args COMPONENTS OPTIONAL_COMPONENTS)
  cmake_parse_arguments(FCLAW "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN} )

  if(NOT DEFINED FCLAW_ROOT)
    set(FCLAW_ROOT ${CMAKE_INSTALL_PREFIX})
  endif()


  # libforestclaw
  if(BUILD_SHARED_LIBS)
    set(FCLAW_LIBRARY_PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
    set(FCLAW_LIBRARY_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(FCLAW_LIBRARY_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
    set(FCLAW_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()


  unset(FCLAW_IMPORTED_LIBARIES)

  add_library(SC::SC STATIC IMPORTED GLOBAL)
  list(APPEND FCLAW_IMPORTED_LIBARIES SC::SC)
  set_target_properties(SC::SC PROPERTIES 
    IMPORTED_LOCATION ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}sc${FCLAW_LIBRARY_SUFFIX})
  target_include_directories(SC::SC INTERFACE ${FCLAW_ROOT}/include)

  add_library(P4EST::P4EST STATIC IMPORTED GLOBAL)
  list(APPEND FCLAW_IMPORTED_LIBARIES P4EST::P4EST)
  set_target_properties(P4EST::P4EST PROPERTIES 
    IMPORTED_LOCATION ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}p4est${FCLAW_LIBRARY_SUFFIX})
  target_include_directories(P4EST::P4EST INTERFACE ${FCLAW_ROOT}/include)
  find_package(ZLIB REQUIRED)
  target_link_libraries(P4EST::P4EST INTERFACE SC::SC ZLIB::ZLIB)

  #FORESTCLAW 
  add_library(FORESTCLAW::FORESTCLAW STATIC IMPORTED GLOBAL)
  list(APPEND FCLAW_IMPORTED_LIBARIES FORESTCLAW::FORESTCLAW)
  set_target_properties(FORESTCLAW::FORESTCLAW PROPERTIES 
    IMPORTED_LOCATION ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}forestclaw${FCLAW_LIBRARY_SUFFIX}
    IMPORTED_LINK_INTERFACE_LANGUAGES  "C;CXX;Fortran")
  target_include_directories(FORESTCLAW::FORESTCLAW INTERFACE ${FCLAW_ROOT}/include)

  target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE P4EST::P4EST)
  if("mpi" IN_LIST FCLAW_COMPONENTS OR "thunderegg" IN_LIST FCLAW_COMPONENTS)
    find_package(MPI COMPONENTS C CXX REQUIRED)
    target_link_libraries(FORESTCLAW::FORESTCLAW INTERFACE MPI::MPI_C INTERFACE MPI::MPI_CXX)
  endif()


  #CLAWPATCH 
  add_library(FORESTCLAW::CLAWPATCH STATIC IMPORTED GLOBAL)
  list(APPEND FCLAW_IMPORTED_LIBARIES FORESTCLAW::CLAWPATCH)
  set_target_properties(FORESTCLAW::CLAWPATCH PROPERTIES 
    IMPORTED_LOCATION ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}clawpatch${FCLAW_LIBRARY_SUFFIX}
    IMPORTED_LINK_INTERFACE_LANGUAGES  "C;CXX;Fortran")
  target_include_directories(FORESTCLAW::CLAWPATCH INTERFACE ${FCLAW_ROOT}/include)

  target_link_libraries(FORESTCLAW::CLAWPATCH INTERFACE FORESTCLAW::FORESTCLAW)


  if("clawpack" IN_LIST FCLAW_COMPONENTS)
    #CLAWPACK4.6
    add_library(FORESTCLAW::CLAWPACK4.6 STATIC IMPORTED GLOBAL)
    list(APPEND FCLAW_IMPORTED_LIBARIES FORESTCLAW::CLAWPACK4.6)
    set_target_properties(FORESTCLAW::CLAWPACK4.6 PROPERTIES 
      IMPORTED_LOCATION ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}clawpack4.6${FCLAW_LIBRARY_SUFFIX}
      IMPORTED_LINK_INTERFACE_LANGUAGES  "C;CXX;Fortran")
    target_include_directories(FORESTCLAW::CLAWPACK4.6 INTERFACE ${FCLAW_ROOT}/include)

    target_link_libraries(FORESTCLAW::CLAWPACK4.6 INTERFACE FORESTCLAW::CLAWPATCH)


    #CLAWPACK5
    add_library(FORESTCLAW::CLAWPACK5 STATIC IMPORTED GLOBAL)
    list(APPEND FCLAW_IMPORTED_LIBARIES FORESTCLAW::CLAWPACK5)
    set_target_properties(FORESTCLAW::CLAWPACK5 PROPERTIES 
      IMPORTED_LOCATION ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}clawpack5${FCLAW_LIBRARY_SUFFIX}
      IMPORTED_LINK_INTERFACE_LANGUAGES  "C;CXX;Fortran")
    target_include_directories(FORESTCLAW::CLAWPACK5 INTERFACE ${FCLAW_ROOT}/include)

    target_link_libraries(FORESTCLAW::CLAWPACK5 INTERFACE FORESTCLAW::CLAWPATCH)


    #CLAWPACK3_46
    add_library(FORESTCLAW::CLAWPACK3_46 STATIC IMPORTED GLOBAL)
    list(APPEND FCLAW_IMPORTED_LIBARIES FORESTCLAW::CLAWPACK3_46)
    set_target_properties(FORESTCLAW::CLAWPACK3_46 PROPERTIES 
      IMPORTED_LOCATION ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}clawpack3_46${FCLAW_LIBRARY_SUFFIX}
      IMPORTED_LINK_INTERFACE_LANGUAGES  "C;CXX;Fortran")
    target_include_directories(FORESTCLAW::CLAWPACK3_46 INTERFACE ${FCLAW_ROOT}/include)

    target_link_libraries(FORESTCLAW::CLAWPACK3_46 INTERFACE FORESTCLAW::CLAWPATCH)
  endif()


  if("thunderegg" IN_LIST FCLAW_COMPONENTS)
    #FC2D_THUNDEREGG
    add_library(FORESTCLAW::FC2D_THUNDEREGG STATIC IMPORTED GLOBAL)
    list(APPEND FCLAW_IMPORTED_LIBARIES FORESTCLAW::FC2D_THUNDEREGG)
    set_target_properties(FORESTCLAW::FC2D_THUNDEREGG PROPERTIES 
      IMPORTED_LOCATION ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}fc2d_thunderegg${FCLAW_LIBRARY_SUFFIX}
      IMPORTED_LINK_INTERFACE_LANGUAGES  "C;CXX;Fortran")
    target_include_directories(FORESTCLAW::FC2D_THUNDEREGG INTERFACE ${FCLAW_ROOT}/include)

    find_package(ThunderEgg CONFIG COMPONENTS P4EST OPTIONAL_COMPONENTS FFTW)
    if(NOT TARGET ThunderEgg::ThunderEgg)
      include(ThunderEggExternalProject)
      ThunderEggExternalProject_GetTarget(COMPONENTS P4EST OPTIONAL_COMPONENTS FFTW)
      list(APPEND FCLAW_IMPORTED_LIBARIES ThunderEgg::ThunderEgg)
      
    endif()

    target_link_libraries(FORESTCLAW::FC2D_THUNDEREGG INTERFACE FORESTCLAW::FORESTCLAW FORESTCLAW::CLAWPATCH ThunderEgg::ThunderEgg)
  endif()


  if("geoclaw" IN_LIST FCLAW_COMPONENTS)
    #GEOCLAW
    add_library(FORESTCLAW::GEOCLAW STATIC IMPORTED GLOBAL)
    list(APPEND FCLAW_IMPORTED_LIBARIES FORESTCLAW::GEOCLAW)
    set_target_properties(FORESTCLAW::GEOCLAW PROPERTIES 
      IMPORTED_LOCATION                  ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}geoclaw${FCLAW_LIBRARY_SUFFIX}
      IMPORTED_LINK_INTERFACE_LANGUAGES  "C;CXX;Fortran")
    target_include_directories(FORESTCLAW::GEOCLAW INTERFACE ${FCLAW_ROOT}/include)

    target_link_libraries(FORESTCLAW::GEOCLAW INTERFACE FORESTCLAW::CLAWPATCH)
  endif()

  if("cudaclaw" IN_LIST FCLAW_COMPONENTS)
    find_package(CUDAToolkit REQUIRED)

    add_library(FORESTCLAW::CUDACLAW STATIC IMPORTED GLOBAL)
    list(APPEND FCLAW_IMPORTED_LIBARIES FORESTCLAW::CUDACLAW)
    set_target_properties(FORESTCLAW::CUDACLAW PROPERTIES 
      IMPORTED_LOCATION                  ${FCLAW_ROOT}/lib/${FCLAW_LIBRARY_PREFIX}cudaclaw${FCLAW_LIBRARY_SUFFIX}
      IMPORTED_LINK_INTERFACE_LANGUAGES  "C;CXX;Fortran;CUDA"
      CUDA_SEPARABLE_COMPILATION ON)
    target_include_directories(FORESTCLAW::CUDACLAW INTERFACE ${FCLAW_ROOT}/include)

    target_link_libraries(FORESTCLAW::CUDACLAW INTERFACE FORESTCLAW::CLAWPATCH CUDA::nvToolsExt)
  endif()

  # --- imported target

  file(MAKE_DIRECTORY ${FCLAW_INCLUDE_DIRS})
  # avoid race condition

  set(FCLAW_IMPORTED_LIBARIES ${FCLAW_IMPORTED_LIBARIES} PARENT_SCOPE)

endfunction(ForestClawExternalProject_GetTarget)

#
# ForestClawExternalProject
# -------------------------
#
# This module defines a function for using ForestClaw as an external project
#
#     ForestClawExternalProject(target
#                               TAG tag
#                               [REPOSITORY repo_url]
#                               [COMPONENTS arg1...]
#                               )
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
#     list of components to build ForestClaw with.
#     Current components are P4EST, FFTW, PETSC, LAPACK
#
function(ForestClawExternalProject)
  include(ExternalProject)
  
  set(options "")
  set(one_value_args REPOSITORY TAG)
  set(multi_value_args COMPONENTS)
  cmake_parse_arguments(FCLAW "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN} )
  string(REPLACE " " ";" FCLAW_COMPONENTS "${FCLAW_COMPONENTS}")
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/local CACHE PATH "install prefix" FORCE)
  endif()
  if(NOT DEFINED FCLAW_ROOT)
    set(FCLAW_ROOT ${CMAKE_INSTALL_PREFIX})
  endif()

  if(NOT DEFINED FCLAW_REPOSITORY)
    set(FCLAW_REPOSITORY "https://github.com/forestclaw/forestclaw")
  endif()

  if(NOT DEFINED FCLAW_TAG)
    message(FATAL_ERROR "TAG for ForestClaw External Project needs to be set")
  endif()

  if(BUILD_SHARED_LIBS)
    set(FCLAW_LIBRARIES ${FCLAW_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ForestClaw${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(FCLAW_LIBRARIES ${FCLAW_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ForestClaw${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()

  set(FCLAW_INCLUDE_DIRS ${FCLAW_ROOT}/include)

  set(FCLAW_CONFIGURE_ARGS -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=${FCLAW_ROOT} -Dapplications=OFF)

  # --- optional dependent libraries
  if("mpi" IN_LIST FCLAW_COMPONENTS OR "thunderegg" IN_LIST FCLAW_COMPONENTS)
    list(APPEND FCLAW_CONFIGURE_ARGS -Dmpi=on)
  endif()

  if("clawpack" IN_LIST FCLAW_COMPONENTS)
    list(APPEND FCLAW_CONFIGURE_ARGS -Dclawpack=on)
  endif()

  if("thunderegg" IN_LIST FCLAW_COMPONENTS)
    list(APPEND FCLAW_CONFIGURE_ARGS -Dthunderegg=on)
  endif()

  if("geoclaw" IN_LIST FCLAW_COMPONENTS)
    list(APPEND FCLAW_CONFIGURE_ARGS -Dgeoclaw=on)
  endif()

  if("cudaclaw" IN_LIST FCLAW_COMPONENTS)
    list(APPEND FCLAW_CONFIGURE_ARGS -Dcudaclaw=on)
  endif()

  list(APPEND FCLAW_CONFIGURE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS})
  list(APPEND FCLAW_CONFIGURE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS})
  list(APPEND FCLAW_CONFIGURE_ARGS -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS})
  list(APPEND FCLAW_CONFIGURE_ARGS -DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS})
  # --- propogate options
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

  # --- imported target

  file(MAKE_DIRECTORY ${FCLAW_INCLUDE_DIRS})
  # avoid race condition

  ForestClawExternalProject_GetTarget(COMPONENTS ${FCLAW_COMPONENTS} OPTIONAL_COMPONENTS ${FCLAW_OPTIONAL_COMPONENTS})

  unset(FCLAW_LIBRARIES)
  foreach(LIBRARY ${FCLAW_IMPORTED_LIBARIES})
    get_target_property(LIBRARY_LOCATION ${LIBRARY} IMPORTED_LOCATION)
    list(APPEND FCLAW_LIBRARIES ${LIBRARY_LOCATION})
  endforeach()

  ExternalProject_Add(ForestClaw
    GIT_REPOSITORY   ${FCLAW_REPOSITORY}
    GIT_TAG          ${FCLAW_TAG}
    CMAKE_ARGS       ${FCLAW_CONFIGURE_ARGS}
    BUILD_BYPRODUCTS ${FCLAW_LIBRARIES}
    DEPENDS          ${FCLAW_EXTERNAL_PROJECT_DEPENDS}
  )

  foreach(LIBRARY ${FCLAW_IMPORTED_LIBARIES})
    add_dependencies(${LIBRARY} ForestClaw)
  endforeach()

endfunction(ForestClawExternalProject)