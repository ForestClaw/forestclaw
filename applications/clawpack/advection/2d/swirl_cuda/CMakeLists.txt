set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(swirl_cuda_f OBJECT 
  psi.f 
  setprob.f
  user_4.6/qinit.f 
  user_4.6/b4step2.f
  user_5.0/qinit.f90 
  user_5.0/b4step2.f90
  user_cuda/qinit.f
  ${all}/clawpack46_setaux.f
  ${all}/clawpack5_setaux.f90
  ${all}/clawpack46_rpn2adv.f 
  ${all}/clawpack46_rpt2adv.f 
  ${all}/clawpack5_rpn2adv.f90 
  ${all}/clawpack5_rpt2adv.f90
)

add_library(swirl_cuda_cu OBJECT 
	user_cuda/rpn2adv.cu
	user_cuda/rpt2adv.cu
	user_cuda/b4step2.cu
)

target_link_libraries(swirl_cuda_cu PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
 FORESTCLAW::CUDACLAW
)

add_executable(swirl_cuda 
  swirl_user.cpp 
  swirl.cpp 
  swirl_options.c
  $<TARGET_OBJECTS:swirl_cuda_f>
  $<TARGET_OBJECTS:swirl_cuda_cu>
)

target_include_directories(swirl_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(swirl_cuda PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
 FORESTCLAW::CUDACLAW
)

# add_test(NAME swirl_cuda
#  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:swirl_cuda>
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
