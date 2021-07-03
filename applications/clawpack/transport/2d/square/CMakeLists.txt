set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(square_f OBJECT 
	value_exceeds_th.f90
	setprob.f
	q0.f
	fdisc.f
	square_basis.f
	mapping_terms.f
	square_velocity.f
	rpn2cons_update_manifold.f
	qexact.f
	dopri5.f
	user_4.6/square_output.f
	user_4.6/compute_error.f
	user_4.6/qinit.f
	user_4.6/square_setaux.f
	user_4.6/rpn2cons_fw_manifold.f
	user_4.6/rpt2cons_manifold.f
	user_4.6/tag4coarsening.f
	user_4.6/tag4refinement.f
	user_5.0/square_output.f
	user_5.0/compute_error.f
	user_5.0/qinit.f
	user_5.0/square_setaux.f
	user_5.0/rpn2cons_fw_manifold.f
	user_5.0/rpt2cons_manifold.f
	user_5.0/tag4coarsening.f
	user_5.0/tag4refinement.f
)

add_executable(square
	square.cpp
	square_user.cpp
	square_options.c
	fclaw2d_map_identity.c
	fclaw2d_map_cart.c
	fclaw2d_map_fivepatch.c
	fclaw2d_map_bilinear.c
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