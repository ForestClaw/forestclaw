set(all ${PROJECT_SOURCE_DIR}/applications/clawpack/advection/3d/all)

add_library(swirl_f OBJECT 
  setprob.f90
  user_4.6/qinit.f90 
  user_4.6/setaux.f90
  ${all}/rpn3_vc_advection.f90 
  ${all}/rpt3_vc_advection.f90 
  ${all}/rptt3_vc_advection.f90 
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
 FORESTCLAW::CLAWPACK3_46
)

add_test(NAME clawpack_advection_3d_swirl COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_3d_swirl PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

