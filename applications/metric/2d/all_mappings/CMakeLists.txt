set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(all_mappings_f OBJECT 
	tag4refinement.f
	metric_output.f
	initialize.f
	compute_error.f
	setprob.f
	diagnostics.f
)

add_executable(all_mappings 
	metric.cpp
	metric_user.cpp
	metric_diagnostics.cpp
	fclaw2d_map_nomap.c
	fclaw2d_map_cart.c
	fclaw2d_map_fivepatch.c
	fclaw2d_map_pillowdisk.c
	fclaw2d_map_pillowdisk5.c
	fclaw2d_map_squareddisk.c
	fclaw2d_map_pillowsphere.c
	fclaw2d_map_cubedsphere.c
	fclaw2d_map_torus.c
  	$<TARGET_OBJECTS:all_mappings_f>
)

target_include_directories(all_mappings PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(all_mappings PRIVATE 
 FORESTCLAW::CLAWPACK4.6
)