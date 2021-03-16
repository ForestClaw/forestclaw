# provides imported target P4EST::P4EST
include(${CMAKE_CURRENT_LIST_DIR}/GitSubmodule.cmake)

set(p4est_external true CACHE BOOL "build p4est library" FORCE)

git_submodule("${PROJECT_SOURCE_DIR}/p4est")
