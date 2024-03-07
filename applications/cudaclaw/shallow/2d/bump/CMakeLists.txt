#set(rp ${PROJECT_SOURCE_DIR}/clawpack/shallow/2d/rp)
set(rpclaw ${PROJECT_SOURCE_DIR}/clawpack/shallow/2d/rp)

add_library(bump_cpu_f OBJECT 
  setprob.f 
  rpn2cons_update.f
  user_4.6/qinit.f90 
  ${rpclaw}/clawpack46_rpn2.f90 
  ${rpclaw}/clawpack46_rpt2.f90 
)

add_library(bump_cuda_f OBJECT 
  user_cuda/rpn2shallow.cu
)

add_executable(bump_cuda 
  bump.cpp 
  bump_options.c 
  bump_user.cpp 
  $<TARGET_OBJECTS:bump_cuda_f>
  $<TARGET_OBJECTS:bump_cpu_f>
)

target_include_directories(bump_cuda PRIVATE 
  FORESTCLAW::CUDACLAW
  FORESTCLAW::CLAWPACK4.6
  FORESTCLAW::CLAWPACK5
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(bump_cuda PRIVATE 
  FORESTCLAW::CUDACLAW
  FORESTCLAW::CLAWPACK4.6
  FORESTCLAW::CLAWPACK5
)

target_include_directories(bump_cuda_f PRIVATE 
  FORESTCLAW::CUDACLAW
  FORESTCLAW::CLAWPACK4.6
  FORESTCLAW::CLAWPACK5
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(bump_cuda_f PRIVATE 
  FORESTCLAW::CUDACLAW
  FORESTCLAW::CLAWPACK4.6
  FORESTCLAW::CLAWPACK5
)


add_test(NAME clawpack_shallow_2d_bump COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_shallow_2d_bump PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

add_test(NAME cudaclaw_shallow_2d_bump COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(cudaclaw_shallow_2d_bump PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")