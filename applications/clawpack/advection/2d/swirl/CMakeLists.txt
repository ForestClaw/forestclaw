set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(swirl_f OBJECT 
  psi.f 
  setprob.f
  user_4.6/qinit.f 
  user_4.6/b4step2.f
  user_5.0/qinit.f90 
  user_5.0/b4step2.f90
  ${all}/clawpack46_setaux.f
  ${all}/clawpack5_setaux.f90
  ${all}/clawpack46_rpn2adv.f 
  ${all}/clawpack46_rpt2adv.f 
  ${all}/clawpack5_rpn2adv.f90 
  ${all}/clawpack5_rpt2adv.f90
)

add_executable(swirl 
  swirl_user.cpp 
  swirl.cpp 
  swirl_options.c
  $<TARGET_OBJECTS:swirl_f>
)

target_include_directories(swirl PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(swirl PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_swirl COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_swirl PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")