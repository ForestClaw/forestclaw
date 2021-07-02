set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(replicated_f OBJECT 
  setprob.f 
  psi.f 
  fdisc.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  ${all}/clawpatch46_tag4refinement.f 
  ${all}/clawpatch46_tag4coarsening.f 
  ${all}/clawpatch5_tag4refinement.f 
  ${all}/clawpatch5_tag4coarsening.f
  ${all}/advection_user_fort.h 
  ${all}/clawpack46_rpn2adv.f 
  ${all}/clawpack46_rpt2adv.f 
  ${all}/clawpack5_rpn2adv.f90 
  ${all}/clawpack5_rpt2adv.f90 
  ${all}/clawpack46_setaux.f 
  ${all}/clawpack5_setaux.f90
)

add_executable(replicated 
  replicated.cpp 
  replicated_options.c 
  replicated_user.cpp 
  $<TARGET_OBJECTS:replicated_f>
)

target_include_directories(replicated PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(replicated PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_replicated COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_replicated PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")