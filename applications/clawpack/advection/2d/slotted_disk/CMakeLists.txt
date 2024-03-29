set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(slotted_disk_f OBJECT 
  fdisc.f 
  slotteddisk.f 
  slotted_disk_exceeds_th.f90
  user_4.6/qinit.f 
  user_5.0/qinit.f90
  ${all}/NCAR_tracer_utils.f 
  ${all}/user46_b4step2_manifold.f 
  ${all}/user5_b4step2_manifold.f90 
  ${all}/user46_setaux_manifold.f 
  ${all}/clawpack46_rpn2adv_manifold.f 
  ${all}/clawpack46_rpt2adv_manifold.f 
  ${all}/user5_setaux_manifold.f90 
  ${all}/clawpack5_rpn2adv_manifold.f90 
  ${all}/clawpack5_rpt2adv_manifold.f90
)

add_executable(slotted_disk 
  slotted_disk.cpp 
  slotted_disk_user.cpp 
  slotted_disk_options.c
  ${all}/advection_patch_setup_manifold.c
  ${all}/advection_b4step2_manifold.c
  ${all}/fclaw2d_map_cubedsphere.c 
  ${all}/fclaw2d_map_pillowsphere.c 
  $<TARGET_OBJECTS:slotted_disk_f>
)

target_include_directories(slotted_disk PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(slotted_disk PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_slotted_disk COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_slotted_disk PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
