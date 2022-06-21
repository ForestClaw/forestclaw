set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(filament_swirl_f OBJECT 
  psi.f 
  filament/psi.f 
  filament/fdisc.f 
  filament/setprob.f 
  filament/user_4.6/qinit.f 
  filament/user_5.0/qinit.f90 
  swirl/psi.f 
  swirl/setprob.f
  swirl/user_4.6/qinit.f 
  swirl/user_4.6/b4step2.f
  swirl/user_5.0/qinit.f90 
  swirl/user_5.0/b4step2.f90
  ${all}/clawpack46_setaux.f 
  ${all}/clawpack5_setaux.f90 
  ${all}/user46_setaux_manifold.f 
  ${all}/user5_setaux_manifold.f90 
  ${all}/clawpack46_rpn2adv.f 
  ${all}/clawpack46_rpt2adv.f 
  ${all}/clawpack46_rpn2adv_manifold.f 
  ${all}/clawpack46_rpt2adv_manifold.f 
  ${all}/clawpack5_rpn2adv.f90 
  ${all}/clawpack5_rpt2adv.f90 
  ${all}/clawpack5_rpn2adv_manifold.f90 
  ${all}/clawpack5_rpt2adv_manifold.f90
)

add_executable(filament_swirl 
  filament_swirl.cpp 
  filament/filament.cpp
  filament/filament_options.c 
  filament/filament_user.cpp 
  swirl/swirl.cpp 
  swirl/swirl_user.cpp 
  swirl/swirl_options.c
  ${all}/advection_patch_setup_manifold.c
  ${all}/fclaw2d_map_cart.c 
  ${all}/fclaw2d_map_fivepatch.c 
  $<TARGET_OBJECTS:filament_swirl_f>
)

target_include_directories(filament_swirl PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(filament_swirl PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_filament_swirl COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_filament_swirl PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")