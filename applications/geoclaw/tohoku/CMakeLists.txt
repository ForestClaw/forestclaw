add_executable(tohoku 
	tohoku_user.cpp
	tohoku_user.h
	tohoku.cpp
)

target_include_directories(tohoku PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(tohoku PRIVATE 
  FORESTCLAW::GEOCLAW
)