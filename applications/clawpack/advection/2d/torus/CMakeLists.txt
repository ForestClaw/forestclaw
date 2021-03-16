set(all ${PROJECT_SOURCE_DIR}/applications/clawpack/advection/2d/all)

add_library(torus_f OBJECT 
  setprob.f 
  torus_velocity.f 
  torus_psi.f 
  fdisc.f 
  user_4.6/qinit.f 
  user_4.6/torus46_setaux.f 
  user_4.6/rpn2cons_fw_manifold.f 
  user_4.6/rpt2cons_manifold.f 
  user_4.6/rpn2cons_update_manifold.f 
  ${all}/clawpack46_rpn2adv_manifold.f 
  ${all}/clawpack46_rpt2adv_manifold.f 
  user_5.0/qinit.f90 
  user_5.0/torus5_setaux.f90 
  ${all}/clawpack5_rpn2adv_manifold.f90 
  ${all}/clawpack5_rpt2adv_manifold.f90
)

add_executable(torus 
  torus.cpp 
  torus_options.c 
  torus_user.cpp 
  fclaw2d_map_torus.c 
  fclaw2d_map_cart.c 
  $<TARGET_OBJECTS:torus_f>
)

target_include_directories(torus PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(torus PRIVATE 
  forestclaw
  clawpack4.6
  clawpack5
  clawpatch
  P4EST::P4EST 
  SC::SC
)