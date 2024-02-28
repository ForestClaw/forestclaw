set(rp ${PROJECT_SOURCE_DIR}/clawpack/shallow/2d/rp)

add_library(radialdam_f OBJECT 
  setprob.f 
  fdisc.f 
  rpn2cons_update_manifold.f
  rpn2cons_update.f
  rpn2cons_update_zero.f
  user_4.6/qinit.f 
  user_4.6/setaux.f90
  user_5.0/qinit.f90 
  user_5.0/setaux.f90 
  ${rp}/clawpack46_rpn2.f90 
  ${rp}/clawpack46_rpt2.f90 
  ${rp}/clawpack46_rpn2_manifold.f90 
  ${rp}/clawpack46_rpt2_manifold.f90 
  ${rp}/clawpack5_rpn2.f90 
  ${rp}/clawpack5_rpt2.f90 
  ${rp}/clawpack5_rpn2_manifold.f90 
  ${rp}/clawpack5_rpt2_manifold.f90)

add_executable(radialdam 
  radialdam.cpp 
  radialdam_options.c 
  radialdam_user.cpp 
  fclaw2d_map_fivepatch.c   
  fclaw2d_map_pillowdisk.c   
  fclaw2d_map_pillowdisk5.c 
  $<TARGET_OBJECTS:radialdam_f>
)

target_include_directories(radialdam PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(radialdam PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_shallow_2d_radialdam COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_shallow_2d_radialdam PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")