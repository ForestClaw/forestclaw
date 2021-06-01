set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(torus_f OBJECT 
	setprob.f
	torus_velocity.f
	torus_psi.f
	dopri5.f
	fdisc.f
	qexact.f
	user_4.6/qinit.f
	user_4.6/compute_error.f
	user_4.6/torus46_setaux.f
	user_4.6/torus46_output.f
	user_4.6/tag4refinement.f
	user_4.6/tag4coarsening.f
	user_4.6/rpn2cons_fw_manifold.f
	user_4.6/rpt2cons_manifold.f
	user_4.6/rpn2cons_update_manifold.f
	${all}/clawpack46_rpn2adv_manifold.f
	${all}/clawpack46_rpt2adv_manifold.f
)

add_executable(torus 
	torus.cpp
	torus_options.c
	torus_user.cpp
	fclaw2d_map_torus.c
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