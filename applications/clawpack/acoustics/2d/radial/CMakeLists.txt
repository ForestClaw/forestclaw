set(rp ${PROJECT_SOURCE_DIR}/clawpack/acoustics/2d/rp)

add_library(radial_f OBJECT 
  setprob.f 
  user_4.6/qinit.f 
  user_4.6/setaux.f90
  user_5.0/qinit.f90 
  user_5.0/setaux.f90 
  ${rp}/clawpack46_rpn2_acoustics.f 
  ${rp}/clawpack46_rpt2_acoustics.f 
  ${rp}/clawpack46_rpn2_acoustics_manifold.f90 
  ${rp}/clawpack46_rpt2_acoustics_manifold.f90 
  ${rp}/clawpack5_rpn2_acoustics.f90 
  ${rp}/clawpack5_rpt2_acoustics.f90 
  ${rp}/clawpack5_rpn2_acoustics_manifold.f90 
  ${rp}/clawpack5_rpt2_acoustics_manifold.f90
)

add_executable(radial 
  radial.cpp 
  radial_user.cpp 
  radial_options.c 
  fclaw2d_map_pillowdisk5.c 
  fclaw2d_map_pillowdisk.c 
  $<TARGET_OBJECTS:radial_f>
)

target_include_directories(radial PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(radial PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_regressions(regressions.txt)
