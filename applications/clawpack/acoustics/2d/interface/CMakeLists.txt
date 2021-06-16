set(rp ${PROJECT_SOURCE_DIR}/clawpack/acoustics/2d/rp)

add_library(interface_f OBJECT 
  p0.f 
  setprob.f 
  value_exceeds_th.f90
  user_4.6/setaux.f 
  user_4.6/qinit.f 
  user_4.6/tag4refinement.f 
  user_4.6/tag4coarsening.f 
  user_5.0/setaux.f90 
  user_5.0/qinit.f90 
  user_5.0/tag4refinement.f90 
  user_5.0/tag4coarsening.f90 
  ${rp}/rpn2acv.f 
  ${rp}/rpt2acv.f 
  ${rp}/rpn2_vc_acoustics.f90 
  ${rp}/rpt2_vc_acoustics.f90
)

add_executable(interface 
  interface.cpp 
  interface_options.c 
  interface_user.cpp 
  $<TARGET_OBJECTS:interface_f>
)

target_include_directories(interface PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(interface PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)