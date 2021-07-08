set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(periodic_f OBJECT 
  setprob.f
  fdisc.f
  user_4.6/qinit.f 
  user_4.6/clawpack46_setaux.f
  user_5.0/qinit.f90 
  user_5.0/clawpack5_setaux.f90
  ${all}/clawpack46_rpn2adv.f 
  ${all}/clawpack46_rpt2adv.f 
  ${all}/clawpack5_rpn2adv.f90 
  ${all}/clawpack5_rpt2adv.f90
)

add_executable(periodic 
  periodic_user.cpp 
  periodic.cpp 
  periodic_options.c
  $<TARGET_OBJECTS:periodic_f>
)

target_include_directories(periodic PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(periodic PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_periodic COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_periodic PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")