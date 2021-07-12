set(all ${PROJECT_SOURCE_DIR}/clawpack/transport/2d/all)

add_library(sphere_transport_f OBJECT 
	setprob.f 
	q0.f 
	fdisc.f 
	sphere_basis.f 
	mapping_terms.f 
	sphere_velocity.f 
	qexact.f 
	dopri5.f 
	user_4.6/sphere46_qinit.f 
	user_4.6/sphere46_output.f 
	user_4.6/sphere46_compute_error.f 
	user_5.0/sphere5_qinit.f 
	user_5.0/sphere5_output.f 
	user_5.0/sphere5_compute_error.f 
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

add_executable(sphere_transport
	sphere.cpp 
	sphere_user.cpp 
	sphere_options.c 
	${all}/transport_patch_setup_manifold.c
	${all}/transport_b4step2_manifold.c
	${all}/fclaw2d_map_cubedsphere.c 
	${all}/fclaw2d_map_pillowsphere.c 
  $<TARGET_OBJECTS:sphere_transport_f>
)

target_include_directories(sphere_transport PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(sphere_transport PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_transport_2d_sphere 
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_transport_2d_sphere
         PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
