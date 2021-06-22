set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(swirlcons_f OBJECT 
	setprob.f
	fdisc.f
	psi.f
	velocity.f
	qexact.f
	dopri5.f
	user_4.6/swirl46_output.f
	user_4.6/compute_error.f
	user_4.6/qinit.f
	user_4.6/swirl46_setaux.f
	user_4.6/rpn2cons_wd.f
	user_4.6/rpn2cons_fw.f
	user_4.6/rpn2cons_ec.f
	user_4.6/rpt2cons.f
	user_4.6/rpn2cons_fw_manifold.f
	user_4.6/rpn2cons_ec_manifold.f
	user_4.6/rpt2cons_manifold.f
	user_4.6/rpn2cons_update_manifold.f
	${all}/clawpack46_rpn2adv_manifold.f
	${all}/clawpack46_rpt2adv_manifold.f 
	user_4.6/tag4coarsening.f
	user_4.6/tag4refinement.f
)

add_executable(swirlcons 
	swirlcons.cpp
	swirlcons_user.cpp
	swirlcons_options.c
	fclaw2d_map_cart.c
	fclaw2d_map_fivepatch.c
	fclaw2d_map_bilinear.c
  $<TARGET_OBJECTS:swirlcons_f>
)

target_include_directories(swirlcons PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(swirlcons PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)