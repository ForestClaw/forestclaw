# set(rp ${PROJECT_SOURCE_DIR}/cudaclaw/euler/2d/rp)

add_library(shockbubble_cuda_f OBJECT 
  fdisc.f 
  setprob.f90 
  user_4.6/qinit.f90
  user_4.6/bc2.f 
  user_4.6/src2.f 
)

add_executable(shockbubble_cuda 
  shockbubble.cpp 
  shockbubble_options.c 
  shockbubble_user.cpp 
	user_cuda/rp2euler.cu
  $<TARGET_OBJECTS:shockbubble_cuda_f>
)

target_include_directories(shockbubble_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(shockbubble_cuda PRIVATE 
 FORESTCLAW::CUDACLAW
)