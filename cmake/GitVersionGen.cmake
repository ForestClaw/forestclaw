# This is an equivalent of build-aux/git-version-gen used in auto tools.

# This CMake script emulates git describe for repositories with a .tarball-version or .repository-version file.
# If the .tarball-version file exists, it reads the version from that file (removing the leading character) and sets the build_number variable accordingly.
# If the .tarball-version file does not exist, it generates a build number based on the version number and commit hash.
# If the current commit is tagged with the version number, the build number will be the version number.
# Otherwise, the build number will be the version number appended with the short commit hash.
# If the repository is dirty, it appends -dirty to the build_number variable.
# This script is intended for use with CMake and shallow copies of repositories.

# Check if the .tarball-version file exists
if(EXISTS "${CMAKE_SOURCE_DIR}/.tarball-version")
  # If .tarball-version exists, read the version from it and remove the leading character
  file(READ "${CMAKE_SOURCE_DIR}/.tarball-version" FORESTCLAW_FULL_VERSION)
  string(STRIP "${FORESTCLAW_FULL_VERSION}" FORESTCLAW_FULL_VERSION)
else()
  # If .tarball-version does not exist, read the version number from the .repository-version file and remove the leading character
  file(READ "${CMAKE_SOURCE_DIR}/.repository-version" FORESTCLAW_FULL_VERSION)
  string(STRIP "${FORESTCLAW_FULL_VERSION}" FORESTCLAW_FULL_VERSION)

  # Get the current commit hash (short version) and check if the current commit is tagged with the version number
  find_package(Git)
  if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD OUTPUT_VARIABLE current_commit OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-list --abbrev-commit -n 1 "refs/tags/${FORESTCLAW_FULL_VERSION}" OUTPUT_VARIABLE tagged_commit OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NOT "${tagged_commit}" STREQUAL "${current_commit}")
      # If the current commit is not tagged, append the short commit hash to the version number
      set(FORESTCLAW_FULL_VERSION "${FORESTCLAW_FULL_VERSION}-${current_commit}")
    endif()

    # Check if the repository is dirty
    execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain -uno -z OUTPUT_VARIABLE git_status)

    if(NOT "${git_status}" STREQUAL "")
        # If the repository is dirty, append -dirty to the build_number variable
        set(FORESTCLAW_FULL_VERSION "${FORESTCLAW_FULL_VERSION}-dirty")
    endif()
  endif()
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

list(GET FORESTCLAW_VERSION_PARTS 0 FORESTCLAW_VERSION_MAJOR)
list(GET FORESTCLAW_VERSION_PARTS 1 FORESTCLAW_VERSION_MINOR)
list(GET FORESTCLAW_VERSION_PARTS 2 FORESTCLAW_VERSION_PATCH)

set(FORESTCLAW_VERSION "${FORESTCLAW_VERSION_MAJOR}.${FORESTCLAW_VERSION_MINOR}.${FORESTCLAW_VERSION_PATCH}")

message (STATUS "ForestClaw version: ${FORESTCLAW_FULL_VERSION}")
