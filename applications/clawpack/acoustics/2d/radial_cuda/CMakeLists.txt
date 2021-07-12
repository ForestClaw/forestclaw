set(rp ${PROJECT_SOURCE_DIR}/clawpack/acoustics/2d/rp)

add_library(radial_cuda_f OBJECT 
  setprob.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  user_5.0/tag4refinement.f90 
  user_5.0/tag4coarsening.f90 
  user_5.0/setaux.f90 
  ${rp}/rpn2ac.f 
	user_cuda/qinit.f
  ${rp}/rpt2ac.f 
  ${rp}/rpn2_acoustics.f90 
  ${rp}/rpt2_acoustics.f90 
  ${rp}/rpn2_acoustics_mapped.f90 
  ${rp}/rpt2_acoustics_mapped.f90
)

add_executable(radial_cuda 
  radial.cpp 
  radial_user.cpp 
  radial_options.c 
  fclaw2d_map_pillowdisk5.c 
	user_cuda/rpn2acoustics.cu
  $<TARGET_OBJECTS:radial_cuda_f>
)

target_include_directories(radial_cuda PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(radial_cuda PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
 FORESTCLAW::CUDACLAW
)