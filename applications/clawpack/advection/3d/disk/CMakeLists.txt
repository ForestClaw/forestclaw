set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/3d/all)
set(all2d ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(disk3dx_f OBJECT 
	psi.f 
	fdisc.f 
	setprob.f 
	user_4.6/qinit.f90 
  ${all}/claw3_setaux_manifold.f90
  ${all}/rpn3_vc_advection.f90 
  ${all}/rpt3_vc_advection.f90 
  ${all}/rptt3_vc_advection.f90 
)

add_executable(disk3dx 
  disk.cpp 
	disk_options.c 
	disk_user.cpp 
  ${all}/advection_patch_setup_manifold.c
	${all2d}/fclaw2d_map_pillowdisk5.c
	${all2d}/fclaw2d_map_pillowdisk.c
  $<TARGET_OBJECTS:disk3dx_f>
)

set_target_properties(disk3dx PROPERTIES OUTPUT_NAME disk)

target_include_directories(disk PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(disk3dx PRIVATE 
 FORESTCLAW::CLAWPACK3_46
)

add_test(NAME clawpack_advection_3d_disk COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_3d_disk PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")


