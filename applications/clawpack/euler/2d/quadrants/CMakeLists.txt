set(rp ${PROJECT_SOURCE_DIR}/clawpack/euler/2d/rp)

add_library(quadrants_f OBJECT 
  setprob.f 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  ${rp}/clawpack46_rpn2_euler4.f 
  ${rp}/clawpack46_rpt2_euler4.f 
  ${rp}/clawpack5_rpn2_euler4.f90 
  ${rp}/clawpack5_rpt2_euler4.f90
)

add_executable(quadrants 
  quadrants.cpp 
  quadrants_options.c 
  quadrants_user.cpp 
  $<TARGET_OBJECTS:quadrants_f>
)

target_include_directories(quadrants PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(quadrants PRIVATE 
  clawpack4.6
  clawpack5
)