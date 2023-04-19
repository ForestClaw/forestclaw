# This is an equivalent of build-aux/git-version-gen used in auto tools.
# It is used to generate the version string for the build.
# It first tries to use git describe.
# If that fails, it uses the version string from the .tarball-version file.
# If all else fails, it sets the version to UNKOWN.

# Try to get the version using git describe
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe  --tags --abbrev=4 --match=v* --dirty=-dirty
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE FORESTCLAW_FULL_VERSION
    RESULT_VARIABLE GIT_RESULT
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if(NOT GIT_RESULT EQUAL "0")
    set(FORESTCLAW_FULL_VERSION "UNKNOWN")
  endif()
else()
  set(FORESTCLAW_FULL_VERSION "UNKNOWN")
endif()


if(EXISTS "${CMAKE_SOURCE_DIR}/.tarball-version" AND FORESTCLAW_FULL_VERSION STREQUAL "UNKNOWN")
  # Read the version from the file
  file(READ "${CMAKE_SOURCE_DIR}/.tarball-version" FORESTCLAW_FULL_VERSION)
  string(STRIP "${FORESTCLAW_FULL_VERSION}" FORESTCLAW_FULL_VERSION)

endif()

# Remove the 'v' prefix from the version string (if present)
string(REGEX REPLACE "^v" "" FORESTCLAW_FULL_VERSION "${FORESTCLAW_FULL_VERSION}")

# Parse the version string into major, minor, and patch variables
string(REGEX MATCHALL "[0-9]+" FORESTCLAW_VERSION_PARTS "${FORESTCLAW_FULL_VERSION}")

# if list size is less than 3, set all verison numbers to 0
list(LENGTH FORESTCLAW_VERSION_PARTS FORESTCLAW_VERSION_PARTS_LENGTH)
if(FORESTCLAW_VERSION_PARTS_LENGTH LESS 3)
  set(FORESTCLAW_VERSION_PARTS 0 0 0)
endif()


set(FORESTCLAW_VERSION "${FORESTCLAW_VERSION_MAJOR}.${FORESTCLAW_VERSION_MINOR}.${FORESTCLAW_VERSION_PATCH}")