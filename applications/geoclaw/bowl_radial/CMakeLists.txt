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