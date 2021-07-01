set(rp ${PROJECT_SOURCE_DIR}/clawpack/shallow/2d/rp)

add_library(bump_f OBJECT 
  setprob.f 
  bump_exceeds_th.f90
  rpn2cons_update.f
  user_4.6/qinit.f 
  user_4.6/tag4refinement.f 
  user_4.6/tag4coarsening.f 
  user_5.0/qinit.f90 
  user_5.0/setaux.f90 
  user_5.0/tag4refinement.f90 
  user_5.0/tag4coarsening.f90 
  ${rp}/clawpack46_rpn2.f 
  ${rp}/clawpack46_rpt2.f 
  ${rp}/clawpack5_rpn2.f90 
  ${rp}/clawpack5_rpt2.f90 
  ${rp}/clawpack5_rpn2_manifold.f90 
  ${rp}/clawpack5_rpt2_manifold.f90
)

add_executable(bump 
  bump.cpp 
  bump_options.c 
  bump_user.cpp 
  $<TARGET_OBJECTS:bump_f>
)

target_include_directories(bump PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(bump PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_shallow_2d_bump COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_shallow_2d_bump PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")