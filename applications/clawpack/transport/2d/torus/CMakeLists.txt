set(all ${PROJECT_SOURCE_DIR}/clawpack/transport/2d/all)

add_library(torus_transport_f OBJECT 
  setprob.f 
  torus_velocity.f 
  torus_psi.f 
  fdisc.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90
  ${all}/transport46_setaux_manifold.f
  ${all}/transport46_set_velocity.f
  ${all}/transport5_setaux_manifold.f
  ${all}/transport5_set_velocity.f
  ${all}/clawpack46_rpn2cons_fw_manifold.f 
  ${all}/clawpack46_rpt2cons_manifold.f 
  ${all}/clawpack5_rpn2cons_fw_manifold.f 
  ${all}/clawpack5_rpt2cons_manifold.f 
  ${all}/rpn2cons_update_manifold.f 
)

# link to clawpatch to get "fclaw2d_metric_terms.i"
target_link_libraries(torus_transport_f PRIVATE
  FORESTCLAW::CLAWPATCH
)

add_executable(torus_transport
  torus.cpp 
  torus_options.c 
  torus_user.cpp 
  ${all}/transport_patch_setup_manifold.c
  ${all}/transport_b4step2_manifold.c
  fclaw2d_map_torus.c 
  $<TARGET_OBJECTS:torus_transport_f>
)

set_target_properties(torus_transport PROPERTIES OUTPUT_NAME torus)

target_include_directories(torus_transport PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(torus_transport PRIVATE
  FORESTCLAW::CLAWPACK4.6
  FORESTCLAW::CLAWPACK5
)

add_regressions(regressions.txt)


