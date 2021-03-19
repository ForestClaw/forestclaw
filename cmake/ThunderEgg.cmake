# provides imported target SC::SC
include(${CMAKE_CURRENT_LIST_DIR}/GitSubmodule.cmake)

set(thunderegg_external true CACHE BOOL "build ThunderEgg library" FORCE)

git_submodule("${PROJECT_SOURCE_DIR}/ThunderEgg")
