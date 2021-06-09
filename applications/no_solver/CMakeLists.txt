
add_library(no_solver_f OBJECT 
	user_4.6/initialize.f
)

add_executable(no_solver 
	no_solver.cpp
	no_solver_user.cpp
	no_solver_options.c
	gem3d_output_mesh.c
  	$<TARGET_OBJECTS:no_solver_f>
)

target_include_directories(no_solver PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(no_solver PRIVATE 
 FORESTCLAW::CLAWPATCH
)