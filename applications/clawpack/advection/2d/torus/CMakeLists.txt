set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

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
  user_5.0/qinit.f90 
  user_5.0/torus5_setaux.f90 
  ${all}/clawpack46_rpn2adv_manifold.f 
  ${all}/clawpack46_rpt2adv_manifold.f 
  ${all}/clawpack5_rpn2adv_manifold.f90 
  ${all}/clawpack5_rpt2adv_manifold.f90
)

# link to clawpatch to get "fclaw2d_metric_terms.i"
target_link_libraries(torus_f PRIVATE
  FORESTCLAW::CLAWPATCH
)

add_executable(torus 
  torus.cpp 
  torus_options.c 
  torus_user.cpp 
  ${all}/fclaw2d_map_torus.c 
  $<TARGET_OBJECTS:torus_f>
)

target_include_directories(torus PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(torus PRIVATE
  FORESTCLAW::CLAWPACK4.6
  FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_torus COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_torus PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")