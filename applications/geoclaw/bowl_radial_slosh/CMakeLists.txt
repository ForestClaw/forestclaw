add_executable(bowl_radial_slosh 
	radial/radial_user.cpp
	radial/radial_options.c
	slosh/slosh_user.cpp
	slosh/slosh_options.c
	bowl.cpp
)

target_include_directories(bowl_radial_slosh PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(bowl_radial_slosh PRIVATE 
  FORESTCLAW::GEOCLAW
)

add_test(NAME geoclaw_bowl_radial_slosh COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoclaw_bowl_radial_slosh PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
