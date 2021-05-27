set(all ${PROJECT_SOURCE_DIR}/applications/clawpack/advection/2d/all)

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
  clawpack4.6
  clawpack5
)

# add_test(NAME swirl
#  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:swirl>
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
