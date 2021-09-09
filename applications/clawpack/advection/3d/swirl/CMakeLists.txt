set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/3d/all)

add_library(swirl3dx_f OBJECT 
  setprob.f90
  user_exceeds_th.f90 
  user_4.6/qinit.f90 
  user_4.6/setaux.f90
  ${all}/rpn3_vc_advection.f90 
  ${all}/rpt3_vc_advection.f90 
  ${all}/rptt3_vc_advection.f90 
)

add_executable(swirl3dx 
  swirl_user.cpp 
  swirl.cpp 
  swirl_options.c
  $<TARGET_OBJECTS:swirl3dx_f>
)
set_target_properties(swirl3dx PROPERTIES OUTPUT_NAME swirl)

target_include_directories(swirl3dx PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(swirl3dx PRIVATE 
 FORESTCLAW::CLAWPACK3_46
)

add_test(NAME clawpack_advection_3d_swirl COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_3d_swirl PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")











