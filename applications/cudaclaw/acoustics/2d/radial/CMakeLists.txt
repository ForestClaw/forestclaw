# set(rp ${PROJECT_SOURCE_DIR}/clawpack/acoustics/2d/rp)

add_library(radial_cuda_f OBJECT 
	user_cuda/qinit.f90
)

add_executable(radial_cuda 
  radial.cpp 
  radial_user.cpp 
  radial_options.c 
	user_cuda/rpn2acoustics.cu
  $<TARGET_OBJECTS:radial_cuda_f>
)

target_include_directories(radial_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(radial_cuda PRIVATE 
 FORESTCLAW::CUDACLAW
)