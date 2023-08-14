set(all2d ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)
set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/3d/all)

add_library(latlong3dx_f OBJECT 
  fdisc.f 
  setprob.f90 
  psi.f90 
  user_4.6/qinit.f90
  ${all}/claw3_setaux_manifold.f90
  ${all}/rpn3_vc_advection.f90 
  ${all}/rpt3_vc_advection.f90 
  ${all}/rptt3_vc_advection.f90 
)

add_executable(latlong_3dx
  latlong.cpp 
  latlong_options.c 
  latlong_user.cpp 
  latlong_map_extrude.c 
  ${all}/advection_patch_setup_manifold.c
  ${all2d}/fclaw2d_map_latlong.c 
  $<TARGET_OBJECTS:latlong3dx_f>
)

set_target_properties(latlong_3dx PROPERTIES OUTPUT_NAME latlong)


target_include_directories(latlong_3dx PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
  ${all2d}
)

target_link_libraries(latlong_3dx PRIVATE 
 FORESTCLAW::CLAWPACK3_46
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)

add_regressions(regressions.txt)
