add_executable(bowl_radial 
	bowl_user.cpp
	bowl.cpp
)

target_include_directories(bowl_radial PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(bowl_radial PRIVATE 
  FORESTCLAW::GEOCLAW
)

add_test(NAME geoclaw_bowl_radial COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(geoclaw_bowl_radial PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
