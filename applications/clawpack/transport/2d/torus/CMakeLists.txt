set(all ${PROJECT_SOURCE_DIR}/clawpack/transport/2d/all)

add_library(torus_transport_f OBJECT 
  setprob.f 
  torus_velocity.f 
  torus_psi.f 
  fdisc.f 
  user_4.6/qinit.f 
  ${all}/clawpatch46_setaux_manifold.f
  ${all}/transport46_set_velocity.f
  ${all}/clawpack46_rpn2cons_fw_manifold.f 
  ${all}/clawpack46_rpt2cons_manifold.f 
  ${all}/rpn2cons_update_manifold.f 
  user_5.0/qinit.f90 
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

target_include_directories(torus_transport PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(torus_transport PRIVATE
  FORESTCLAW::CLAWPACK4.6
  FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_transport_2d_torus COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_transport_2d_torus PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

