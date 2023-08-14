set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/3d/all)
set(all2d ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)


add_library(swirl3dx_f OBJECT 
  setprob.f90
  fdisc.f
  user_exceeds_th.f90 
  user_4.6/qinit.f90 
  user_4.6/swirl_setaux.f90
  user_4.6/swirl_setaux_manifold.f90
  ${all}/rpn3_vc_advection.f90 
  ${all}/rpt3_vc_advection.f90 
  ${all}/rptt3_vc_advection.f90 
)

add_executable(swirl3dx 
  swirl_user.cpp 
  swirl.cpp 
  swirl_options.c
  swirl_map_extrude.c
  ${all2d}/fclaw2d_map_cart.c 
  ${all2d}/fclaw2d_map_fivepatch.c 
  ${all2d}/fclaw2d_map_bilinear.c 
  $<TARGET_OBJECTS:swirl3dx_f>
)

set_target_properties(swirl3dx PROPERTIES OUTPUT_NAME swirl)

target_include_directories(swirl3dx PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(swirl3dx PRIVATE 
 FORESTCLAW::CLAWPACK3_46
)

add_regressions(regressions.txt)










