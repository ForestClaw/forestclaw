set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(correlatedcb_f OBJECT 
  cosine_bell.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  ${all}/NCAR_tracer_utils.f 
  ${all}/user46_b4step2_manifold.f 
  ${all}/user5_b4step2_manifold.f90 
  ${all}/user46_setaux_manifold.f 
  ${all}/user5_setaux_manifold.f90 
  ${all}/clawpack46_rpn2adv_manifold.f 
  ${all}/clawpack46_rpt2adv_manifold.f 
  ${all}/clawpack5_rpn2adv_manifold.f90 
  ${all}/clawpack5_rpt2adv_manifold.f90
)

add_executable(correlatedcb 
  correlatedcb.cpp 
  correlatedcb_user.cpp 
  correlatedcb_options.c 
  ${all}/advection_patch_setup_manifold.c
  ${all}/advection_b4step2_manifold.c
  ${all}/fclaw2d_map_cubedsphere.c 
  ${all}/fclaw2d_map_pillowsphere.c 
  $<TARGET_OBJECTS:correlatedcb_f>
)

target_include_directories(correlatedcb PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(correlatedcb PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_correlatedcb COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_correlatedcb PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")