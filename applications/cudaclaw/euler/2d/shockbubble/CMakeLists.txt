# set(rp ${PROJECT_SOURCE_DIR}/cudaclaw/euler/2d/rp)
set(rpclaw ${PROJECT_SOURCE_DIR}/clawpack/euler/2d/rp)

add_library(shockbubble_cpu_f OBJECT 
  user_4.6/qinit.f 
  user_4.6/setaux.f 
  user_4.6/bc2.f 
  user_4.6/src2.f 
  ${rpclaw}/clawpack46_rpn2_euler4.f 
  ${rpclaw}/clawpack46_rpt2_euler4.f 
  ${rpclaw}/clawpack46_rpn2_euler5.f 
  ${rpclaw}/clawpack46_rpt2_euler5.f 
)

add_library(shockbubble_cuda_f OBJECT 
  fdisc.f 
  setprob.f90 
  user_cuda/qinit.f90
  user_cuda/bc2.f 
  user_cuda/src2.f 
)

add_executable(shockbubble_cuda 
  shockbubble.cpp 
  shockbubble_options.c 
  shockbubble_user.cpp 
	user_cuda/rp2euler.cu
  $<TARGET_OBJECTS:shockbubble_cuda_f>
  $<TARGET_OBJECTS:shockbubble_cpu_f>
)

target_include_directories(shockbubble_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(shockbubble_cuda PRIVATE 
 FORESTCLAW::CUDACLAW
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_euler_2d_shockbubble COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_euler_2d_shockbubble PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")