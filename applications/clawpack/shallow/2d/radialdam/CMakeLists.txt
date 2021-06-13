set(rp ${PROJECT_SOURCE_DIR}/clawpack/shallow/2d/rp)

add_library(radialdam_f OBJECT 
  setprob.f 
  fdisc.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  user_5.0/setaux.f90 
  user_5.0/tag4refinement.f90 
  user_5.0/tag4coarsening.f90 
  ${rp}/clawpack46_rpn2.f 
  ${rp}/clawpack46_rpt2.f 
  ${rp}/clawpack5_rpn2.f90 
  ${rp}/clawpack5_rpt2.f90 
  ${rp}/clawpack5_rpn2_manifold.f90 
  ${rp}/clawpack5_rpt2_manifold.f90
)

add_executable(radialdam 
  radialdam.cpp 
  radialdam_options.c 
  radialdam_user.cpp 
  fclaw2d_map_pillowdisk.c 
  $<TARGET_OBJECTS:radialdam_f>
)

target_include_directories(radialdam PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(radialdam PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)