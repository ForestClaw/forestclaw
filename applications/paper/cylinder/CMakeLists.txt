set(all ${PROJECT_SOURCE_DIR}/paper/all)

add_library(cylinder_f OBJECT 
	setprob.f
	compute_exact.f90
	cylinder_compute_area_fort.f90 
	latlong_compute_area_fort.f90
	dopri5.f
	qexact.f
	q0.f
	fdisc.f
	user_4.6/qinit.f
	user_4.6/compute_error.f
	user_4.6/cylinder_setaux.f
	user_4.6/cylinder46_output.f
	user_4.6/cylinder_tag4coarsening.f
	user_4.6/cylinder_tag4refinement.f
	user_4.6/rpn2cons_fw_manifold.f
  user_4.6/rpt2cons_manifold.f
	user_4.6/rpn2cons_update_manifold.f
	${all}/mapc2m_cylinder.f90
	${all}/cylinder_velocity.f90
	${all}/cylinder_basis.f90
	${all}/mapc2m_latlong2.f90
	${all}/latlong_velocity.f90
	${all}/latlong_basis.f90
	${all}/mapping_terms.f
)

add_executable(cylinder
	cylinder.cpp 
	cylinder_user.cpp 
	cylinder_options.c 
	fclaw2d_map_cylinder.c
  $<TARGET_OBJECTS:cylinder_f>
)

target_include_directories(cylinder PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(cylinder PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME paper_cylinder 
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(paper_cylinder
         PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
