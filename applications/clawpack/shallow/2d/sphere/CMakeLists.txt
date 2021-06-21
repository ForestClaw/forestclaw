add_library(sphere_f OBJECT
	setprob.f
	write_latlong_data.f
	user_4.6/qinit.f
	user_4.6/gaussian.f
	user_4.6/setaux_sphere.f
	user_4.6/b4step2.f
	user_4.6/rpn2swsphere.f
	user_4.6/rpt2swq.f
)

add_executable(sphere
	sphere.cpp
	sphere_user.cpp
	fclaw2d_map_latlong.c
	fclaw2d_map_cubedsphere.c
	$<TARGET_OBJECTS:sphere_cuda_f>
)

target_include_directories(sphere PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(sphere PRIVATE 
	FORESTCLAW::CLAWPACK4.6
	FORESTCLAW::CLAWPACK5
)