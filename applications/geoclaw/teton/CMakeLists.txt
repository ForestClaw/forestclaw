add_library(teton_f OBJECT
	reservoir_module.f90
	qinit.f90
)

target_link_libraries(teton_f PRIVATE 
  geoclaw
)

add_executable(teton 
	teton_user.cpp
	teton_user.h
	teton.cpp
)

target_include_directories(teton PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(teton PRIVATE 
  teton_f
  geoclaw
)