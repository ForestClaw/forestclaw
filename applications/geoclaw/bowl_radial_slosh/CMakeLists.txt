add_executable(bowl_radial_slosh 
	radial/radial_user.cpp
	radial/radial_options.c
	slosh/slosh_user.cpp
	slosh/slosh_options.c
	slosh/qinit.f90
	bowl.cpp
)

target_include_directories(bowl_radial_slosh PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(bowl_radial_slosh PRIVATE 
  FORESTCLAW::GEOCLAW
)

add_regressions(regressions.txt)

