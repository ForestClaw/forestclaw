set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(disk_f OBJECT 
	psi.f 
	fdisc.f 
	setprob.f 
	user_4.6/qinit.f 
	user_5.0/qinit.f90 
	${all}/user46_setaux_manifold.f 
	${all}/clawpack46_rpn2adv_manifold.f 
	${all}/clawpack46_rpt2adv_manifold.f 
	${all}/user5_setaux_manifold.f90 
	${all}/clawpack5_rpn2adv_manifold.f90 
	${all}/clawpack5_rpt2adv_manifold.f90 
)

add_executable(disk 
  disk.cpp 
	disk_options.c 
	disk_user.cpp 
	${all}/fclaw2d_map_pillowdisk5.c
	${all}/fclaw2d_map_pillowdisk.c
  $<TARGET_OBJECTS:disk_f>
)

target_include_directories(disk PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(disk PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_test(NAME clawpack_advection_2d_disk COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_2d_disk PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")