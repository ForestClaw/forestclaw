set(rp ${PROJECT_SOURCE_DIR}/applications/clawpack/shallow/2d/rp)

add_library(bump_f OBJECT 
  setprob.f 
  user_4.6/qinit.f 
  user_4.6/tag4refinement.f 
  user_4.6/tag4coarsening.f 
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

add_executable(bump 
  bump.cpp 
  bump_options.c 
  bump_user.cpp 
  fclaw2d_map_fivepatch.c 
  $<TARGET_OBJECTS:bump_f>
)

target_include_directories(bump PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(bump PRIVATE 
  clawpack4.6
  clawpack5
)