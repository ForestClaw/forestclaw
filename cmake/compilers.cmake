if(MSVC)
  add_compile_options(
  "$<$<COMPILE_LANGUAGE:C,CXX>:/W4>"
  )
else()
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:C,CXX>:-Wall>"
  )
endif()

if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:Fortran>:-warn>"
  )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:Fortran>:-Wall>"
    )
endif()

# preprocess all Fortran files
set(CMAKE_Fortran_PREPROCESS ON)