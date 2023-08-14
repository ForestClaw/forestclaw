set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(filament_f OBJECT 
  psi.f 
  fdisc.f 
  setprob.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
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

add_executable(filament 
  filament.cpp 
  filament_options.c 
  filament_user.cpp 
  ${all}/advection_patch_setup_manifold.c
  ${all}/fclaw2d_map_cart.c 
  ${all}/fclaw2d_map_fivepatch.c 
  ${all}/fclaw2d_map_bilinear.c 
  $<TARGET_OBJECTS:filament_f>
)

target_include_directories(filament PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(filament PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_regressions(regressions.txt)
