set(all ${PROJECT_SOURCE_DIR}/applications/clawpack/transport/all)
set(advection_all ${PROJECT_SOURCE_DIR}/applications/clawpack/advection/2d/all)

add_library(slotted_disk_f OBJECT 
  setprob.f 
  fdisc.f 
  slotteddisk.f 
  user_4.6/qinit.f 
  user_4.6/tag4refinement.f 
  user_4.6/tag4coarsening.f 
  user_5.0/qinit.f90 
  user_5.0/tag4refinement.f 
  user_5.0/tag4coarsening.f 
  ${all}/velocity_wind.f 
  ${all}/gaussian.f 
  ${all}/user46_b4step2_manifold.f 
  ${all}/user5_b4step2_manifold.f90 
  ${advection_all}/advection_user_fort.h 
  ${advection_all}/user46_setaux_manifold.f 
  ${advection_all}/clawpack46_rpn2adv_manifold.f 
  ${advection_all}/clawpack46_rpt2adv_manifold.f 
  ${advection_all}/user5_setaux_manifold.f90 
  ${advection_all}/clawpack5_rpn2adv_manifold.f90 
  ${advection_all}/clawpack5_rpt2adv_manifold.f90
)

add_executable(slotted_disk 
  slotted_disk.cpp 
  slotted_disk_user.cpp 
  ${all}/transport_user.cpp 
  ${all}/transport_domain.c 
  ${all}/transport_options.c 
  ${all}/fclaw2d_map_cubedsphere.c 
  ${all}/fclaw2d_map_pillowsphere.c 
  $<TARGET_OBJECTS:slotted_disk_f>
)

target_include_directories(slotted_disk PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
  ${advection_all}
)

target_link_libraries(slotted_disk PRIVATE 
  clawpack4.6
  clawpack5
)