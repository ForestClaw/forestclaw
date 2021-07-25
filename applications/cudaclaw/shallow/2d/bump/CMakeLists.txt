# set(rp ${PROJECT_SOURCE_DIR}/clawpack/shallow/2d/rp)

add_library(bump_cuda_f OBJECT 
  user_cuda/qinit.f 
)

add_executable(bump_cuda 
  bump.cpp 
  bump_options.c 
  bump_user.cpp 
	user_cuda/rpn2shallow.cu
  $<TARGET_OBJECTS:bump_cuda_f>
)

target_include_directories(bump_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(bump_cuda PRIVATE 
 FORESTCLAW::CUDACLAW
)