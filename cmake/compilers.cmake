include(CheckCCompilerFlag)

# --- compiler options

check_c_compiler_flag(-Wall _has_wall)
if(_has_wall)
  add_compile_options(
    $<$<COMPILE_LANGUAGE:C>:-Wall>
    $<$<COMPILE_LANGUAGE:CXX>:-Wall>
    $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
    )
else()
  check_c_compiler_flag(/Wall _has_msvc_wall)
  if(_has_msvc_wall)
    add_compile_options(/Wall)
  endif()
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  # something
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  add_compile_options(
#    $<$<COMPILE_LANGUAGE:Fortran>:-std=f2003>
    $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument>
    $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-label>
    )
endif()