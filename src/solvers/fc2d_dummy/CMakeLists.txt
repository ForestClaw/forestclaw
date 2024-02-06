# -- add library
add_library(dummy)
set_target_properties(dummy PROPERTIES EXPORT_NAME CLAWPACK5)

target_sources(dummy PRIVATE
	fc2d_dummy.cpp
)

target_include_directories(dummy
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

# -- install
install(FILES
	fc2d_dummy.h
  TYPE INCLUDE
)

install(TARGETS dummy
  EXPORT ${PROJECT_NAME}-targets)

# imported target, for use from FetchContent
add_library(FORESTCLAW::DUMMY INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::DUMMY INTERFACE dummy)
