add_executable(torthem 
  torthem.c 
  torus_common.h 
  torus_common.c
)

target_include_directories(torthem PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(torthem PRIVATE 
  forestclaw
  clawpack4.6
  clawpack5
  clawpatch
  P4EST::P4EST 
  SC::SC
)