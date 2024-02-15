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
    "$<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument,-Wno-unused-variable,-Wno-unused-label>"
    )
endif()
