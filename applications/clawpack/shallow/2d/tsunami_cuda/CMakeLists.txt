
set(rp ${PROJECT_SOURCE_DIR}/clawpack/shallow/2d/rp)

add_library(tsunami_cuda_f OBJECT
  setprob.f
  bathy.f90
  user_4.6/qinit.f90
  user_4.6/tag4refinement.f
  user_4.6/tag4coarsening.f
  user_4.6/setaux.f90
  ${rp}/rpn2_tsunami.f90
  ${rp}/rpt2_tsunami.f90
  ${rp}/geoclaw_riemann_utils.f
)

add_executable(tsunami_cuda
  tsunami.cpp
  tsunami_options.c
  tsunami_user.cpp
  user_cuda/rpn2_tsunami.cu
  $<TARGET_OBJECTS:tsunami_cuda_f>
)

target_include_directories(tsunami_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(tsunami_cuda PRIVATE 
	FORESTCLAW::CLAWPACK4.6
	FORESTCLAW::CLAWPACK5
	FORESTCLAW::CUDACLAW
)