set(all ${PROJECT_SOURCE_DIR}/clawpack/transport/2d/all)

add_library(square_f OBJECT 
	setprob.f
	q0.f
	fdisc.f
	square_basis.f
	mapping_terms.f
	square_velocity.f
	qexact.f
	user_4.6/square46_qinit.f
	user_4.6/square46_output.f
	user_4.6/square46_compute_error.f
	user_5.0/square5_qinit.f
	user_5.0/square5_output.f
	user_5.0/square5_compute_error.f
	${all}/transport46_set_velocity.f
	${all}/transport46_setaux_manifold.f
  ${all}/clawpack46_rpn2cons_fw_manifold.f 
  ${all}/clawpack46_rpt2cons_manifold.f 
	${all}/transport5_set_velocity.f
	${all}/transport5_setaux_manifold.f
  ${all}/clawpack5_rpn2cons_fw_manifold.f 
  ${all}/clawpack5_rpt2cons_manifold.f 
  ${all}/rpn2cons_update_manifold.f 
)

add_executable(square
	square.cpp
	square_user.cpp
	square_options.c
	${all}/transport_patch_setup_manifold.c
	${all}/transport_b4step2_manifold.c
	${all}/fclaw2d_map_identity.c
	${all}/fclaw2d_map_cart.c
	${all}/fclaw2d_map_fivepatch.c
	${all}/fclaw2d_map_bilinear.c
  $<TARGET_OBJECTS:square_f>
)

target_include_directories(square PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(square PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_regressions(regressions.txt)

