set(all ${PROJECT_SOURCE_DIR}/demo/2d/all)

add_library(swirl_restart_f OBJECT 
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

add_executable(swirl_restart 
  swirl_user.cpp 
  swirl.cpp 
  swirl_options.c
  swirl_ray.c
  $<TARGET_OBJECTS:swirl_f>
)

set_target_properties(swirl_restart PROPERTIES OUTPUT_NAME swirl)

target_include_directories(swirl_restart PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(swirl_restart PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

