function(add_regressions TARGET)
  if(NOT _WORKING_DIRECTORY)
    set(_WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
  if(NOT _TEST_LIST)
    set(_TEST_LIST ${TARGET}_TESTS)
  endif()

  ## Generate a unique name based on the extra arguments
  string(SHA1 args_hash "${TARGET} ${CMAKE_CURRENT_BINARY_DIR}")
  string(SUBSTRING ${args_hash} 0 7 args_hash)

  # Define rule to generate test list for aforementioned test executable
  set(ctest_include_file "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}_include-${args_hash}.cmake")
  file(WRITE "${ctest_include_file}"
    "include(\"${PROJECT_SOURCE_DIR}/cmake/Modules/discover_regressions.cmake\")\n"
    "add_regressions(\"${CMAKE_CURRENT_SOURCE_DIR}/${TARGET}\" \"${CMAKE_CURRENT_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}\" \"${MPIEXEC_EXECUTABLE}\" \"${MPIEXEC_NUMPROC_FLAG}\" \"${MPIEXEC_MAX_NUMPROCS}\")\n"
  )

  if(NOT CMAKE_VERSION VERSION_LESS 3.10)
    # Add discovered tests to directory TEST_INCLUDE_FILES
    set_property(DIRECTORY
      APPEND PROPERTY TEST_INCLUDE_FILES "${ctest_include_file}"
    )
  else()
    # Add discovered tests as directory TEST_INCLUDE_FILE if possible
    get_property(test_include_file_set DIRECTORY PROPERTY TEST_INCLUDE_FILE SET)
    if(NOT ${test_include_file_set})
      set_property(DIRECTORY
        PROPERTY TEST_INCLUDE_FILE "${ctest_include_file}"
      )
    else()
      message(FATAL_ERROR
        "Cannot set more than one TEST_INCLUDE_FILE"
      )
    endif()
  endif()

endfunction()

###############################################################################

set(_DISCOVER_REGRESSIONS_SCRIPT
  ${CMAKE_CURRENT_LIST_DIR}/discover_regressions.cmake
)
