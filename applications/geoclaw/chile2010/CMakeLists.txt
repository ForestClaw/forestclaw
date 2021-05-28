add_executable(chile2010 
	chile2010.cpp
)

target_include_directories(chile2010 PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(chile2010 PRIVATE 
  FORESTCLAW::GEOCLAW
)