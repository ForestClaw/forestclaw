find_package(Doxygen)
option(build_documentation "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(build_documentation)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()
  if(NOT DEFINED THUNDEREGG_HTML_OUTPUT_DIR)
    set(THUNDEREGG_HTML_OUTPUT_DIR "html")
  endif()

  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  #include(FetchContent)
  #FetchContent_Declare(
  #  doxygencss
  #  GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css
  #  GIT_TAG        main
  #)
  #FetchContent_MakeAvailable(doxygencss)
  configure_file(${doxyfile_in} ${doxyfile} @ONLY)
  add_custom_target(
    doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
endif()