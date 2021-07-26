set(all ${PROJECT_SOURCE_DIR}/cudaclaw/advection/2d/all)

add_library(swirl_cuda_f OBJECT 
  psi.f 
  user_cuda/qinit.f
)

add_executable(swirl_cuda
  swirl_user.cpp 
  swirl.cpp 
  swirl_options.c
	user_cuda/rpn2adv.cu
	user_cuda/rpt2adv.cu
	user_cuda/b4step2.cu
  $<TARGET_OBJECTS:swirl_cuda_f>
)

target_include_directories(swirl_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(swirl_cuda PRIVATE 
 FORESTCLAW::CUDACLAW
)

add_test(NAME cudaclaw_advection_2d_swirl COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(cudaclaw_advection_2d_swirl PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
