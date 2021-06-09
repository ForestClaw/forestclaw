set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(periodic_f OBJECT 
	setprob.f
	${all}/advection_user_fort.h
	user_4.6/qinit.f
	user_4.6/periodic_output.f
	user_4.6/compute_error.f
	user_4.6/limiter3.f
	user_4.6/rpn2.f
	user_4.6/rpt2.f
	user_4.6/rpn2cons_update.f
	user_4.6/flux2.f
	user_4.6/tag4refinement.f
	user_4.6/tag4coarsening.f
	user_4.6/periodic_interpolate.f
	${all}/clawpack5_rpn2adv.f90
	${all}/clawpack5_rpt2adv.f90
)

add_executable(periodic 
	periodic_user.cpp
	periodic_options.c
	periodic.cpp
  $<TARGET_OBJECTS:periodic_f>
)

target_include_directories(periodic PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(periodic PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)