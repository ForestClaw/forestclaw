set(rp ${PROJECT_SOURCE_DIR}/clawpack/acoustics/2d/rp)

add_library(interface_f OBJECT 
  p0.f 
  setprob.f 
  user_4.6/setaux.f 
  user_4.6/qinit.f 
  user_5.0/setaux.f90 
  user_5.0/qinit.f90 
  ${rp}/clawpack46_rpn2_acoustics_vc.f 
  ${rp}/clawpack46_rpt2_acoustics_vc.f 
  ${rp}/clawpack5_rpn2_acoustics_vc.f90 
  ${rp}/clawpack5_rpt2_acoustics_vc.f90
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

add_regressions(regressions.txt)
