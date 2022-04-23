# set(rp ${PROJECT_SOURCE_DIR}/clawpack/acoustics/2d/rp)
set(rpclaw ${PROJECT_SOURCE_DIR}/clawpack/acoustics/2d/rp)

add_library(radial_cpu_f OBJECT 
  setprob.f 
  user_4.6/qinit.f 
  user_4.6/setaux.f90
  ${rpclaw}/clawpack46_rpn2_acoustics.f 
  ${rpclaw}/clawpack46_rpt2_acoustics.f 
)

add_library(radial_cuda_f OBJECT 
	user_cuda/qinit.f90
)

add_executable(radial_cuda 
  radial.cpp 
  radial_user.cpp 
  radial_options.c 
	user_cuda/rpn2acoustics.cu
  $<TARGET_OBJECTS:radial_cuda_f>
  $<TARGET_OBJECTS:radial_cpu_f>
)

target_include_directories(radial_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
  ${rpclaw}
)

target_link_libraries(radial_cuda PRIVATE 
 FORESTCLAW::CUDACLAW
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_acoustics_2d_radial COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_acoustics_2d_radial PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")