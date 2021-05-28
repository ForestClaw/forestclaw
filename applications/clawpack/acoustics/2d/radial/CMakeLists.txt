set(rp ${PROJECT_SOURCE_DIR}/clawpack/acoustics/2d/rp)

add_library(radial_f OBJECT 
  setprob.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  user_5.0/tag4refinement.f90 
  user_5.0/tag4coarsening.f90 
  user_5.0/setaux.f90 
  ${rp}/rpn2ac.f 
  ${rp}/rpt2ac.f 
  ${rp}/rpn2_acoustics.f90 
  ${rp}/rpt2_acoustics.f90 
  ${rp}/rpn2_acoustics_mapped.f90 
  ${rp}/rpt2_acoustics_mapped.f90
)

add_executable(radial 
  radial.cpp 
  radial_user.cpp 
  radial_options.c 
  fclaw2d_map_pillowdisk5.c 
  $<TARGET_OBJECTS:radial_f>
)

target_include_directories(radial PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(radial PRIVATE 
  clawpack4.6
  clawpack5
)