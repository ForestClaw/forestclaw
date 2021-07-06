set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(mesh_f OBJECT 
	psi.f
	setprob.f
	${all}/clawpack46_setaux.f
	${all}/clawpack_user.h
	${all}/clawpack5_setaux.f90
	${all}/clawpack46_rpn2adv.f
	${all}/clawpack46_rpt2adv.f
	${all}/clawpack5_rpn2adv.f90
	${all}/clawpack5_rpt2adv.f90
	user_4.6/qinit.f
	user_4.6/b4step2.f
	user_5.0/qinit.f90
	user_5.0/b4step2.f90
)

add_executable(mesh 
	mesh_user.cpp
	mesh.cpp
	fclaw2d_map_nomap.c
  	$<TARGET_OBJECTS:mesh_f>
)

target_include_directories(mesh PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(mesh PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)