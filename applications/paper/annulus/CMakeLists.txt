add_library(annulus_f OBJECT 
	setprob.f
	fdisc.f
	annulus_psi.f
	annulus_velocity.f
	user_4.6/qinit.f
	user_4.6/b4step2.f90
	user_4.6/rpn2cons_fw_manifold.f
	user_4.6/rpt2cons_manifold.f
	user_4.6/rpn2qad_flux.f90
	user_4.6/annulus46_setaux.f
	user_4.6/tag4refinement.f
	user_4.6/tag4coarsening.f
	user_4.6/annulus46_rpt2adv_manifold.f
)

add_executable(annulus 
	annulus.cpp
	annulus_options.c
	annulus_user.cpp
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