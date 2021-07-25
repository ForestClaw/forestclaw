set(all ${PROJECT_SOURCE_DIR}/paper/all)

add_library(sphere_paper_f OBJECT 
	setprob.f 
	q0.f 
	fdisc.f 
	sphere_basis.f 
	mapping_terms.f 
	sphere_velocity.f 
	qexact.f 
	dopri5.f 
	user_4.6/sphere_output.f
	user_4.6/compute_error.f
	user_4.6/qinit.f
	user_4.6/sphere_setaux.f
	user_4.6/rpn2cons_fw_manifold.f
	user_4.6/rpt2cons_manifold.f
	user_4.6/rpn2_cons_update_manifold.f
	user_4.6/tag4coarsening.f
	user_4.6/tag4refinement.f
)

add_executable(sphere_paper
	sphere.cpp 
	sphere_user.cpp 
	sphere_options.c 
	fclaw2d_map_cubedsphere.c 
	fclaw2d_map_pillowsphere.c 
  $<TARGET_OBJECTS:sphere_paper_f>
)

target_include_directories(sphere_paper PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(sphere_paper PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME paper_sphere 
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(paper_sphere
         PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
