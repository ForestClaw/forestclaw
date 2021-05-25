set(rp ${PROJECT_SOURCE_DIR}/applications/clawpack/burgers/2d/rp)

add_library(pwconst_f OBJECT 
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  ${rp}/clawpack46_rpn2_burgers.f  
  ${rp}/clawpack46_rpt2_burgers.f  
  ${rp}/clawpack5_rpn2_burgers.f90 
  ${rp}/clawpack5_rpt2_burgers.f90
)

add_executable(pwconst 
  pwconst.cpp 
  pwconst_options.c 
  pwconst_user.cpp 
  $<TARGET_OBJECTS:pwconst_f>
)

target_include_directories(pwconst PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(pwconst PRIVATE 
  forestclaw
  clawpack4.6
  clawpack5
  clawpatch
)