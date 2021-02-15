set(a ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(swirl_f OBJECT psi.f setprob.f
user_4.6/qinit.f user_4.6/b4step2.f
user_5.0/qinit.f90 user_5.0/b4step2.f90
${a}/clawpack46_setaux.f ${a}/clawpack5_setaux.f90 ${a}/clawpack46_rpn2adv.f ${a}/clawpack46_rpt2adv.f ${a}/clawpack5_rpn2adv.f90 ${a}/clawpack5_rpt2adv.f90)

add_executable(swirl swirl_user.cpp swirl.cpp swirl_options.c
  $<TARGET_OBJECTS:swirl_f>)
target_include_directories(swirl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${a})
target_link_libraries(swirl PRIVATE FORESTCLAW::FORESTCLAW
FORESTCLAW::clawpack4.6 FORESTCLAW::clawpack5 FORESTCLAW::clawpatch)

add_test(NAME swirl
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:swirl>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
