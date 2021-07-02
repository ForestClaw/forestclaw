set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(annulus_f OBJECT 
  psi.f 
  setprob.f
  fdisc.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  ${all}/user46_setaux_manifold.f
  ${all}/user5_setaux_manifold.f90
  ${all}/clawpack46_rpn2adv_manifold.f 
  ${all}/clawpack46_rpt2adv_manifold.f 
  ${all}/clawpack5_rpn2adv_manifold.f90 
  ${all}/clawpack5_rpt2adv_manifold.f90
  ${all}/clawpatch46_tag4refinement.f 
  ${all}/clawpatch46_tag4coarsening.f 
  ${all}/clawpatch5_tag4refinement.f 
  ${all}/clawpatch5_tag4coarsening.f
)

add_executable(annulus 
  annulus_user.cpp 
  annulus.cpp 
  annulus_options.c
	fclaw2d_map_annulus.c
  $<TARGET_OBJECTS:annulus_f>
)

target_include_directories(annulus PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(annulus PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_annulus COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_annulus PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")