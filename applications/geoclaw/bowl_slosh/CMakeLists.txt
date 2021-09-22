add_library(bowl_slosh_f OBJECT
	qinit.f90
)

target_link_libraries(bowl_slosh_f PRIVATE 
  FORESTCLAW::GEOCLAW
)

add_executable(bowl_slosh 
	slosh_user.cpp
	slosh_user.h
	slosh.cpp
	slosh_options.c
	$<TARGET_OBJECTS:bowl_slosh_f>	
)

target_link_libraries(bowl_slosh_f PRIVATE 
  FORESTCLAW::GEOCLAW
)

target_include_directories(bowl_slosh PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(bowl_slosh PRIVATE 
  FORESTCLAW::GEOCLAW
  bowl_slosh_f
)

add_test(NAME geoclaw_bowl_slosh COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoclaw_bowl_slosh PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")

