set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/3d/all)
set(all2d ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(filament3dx_f OBJECT 
  psi.f 
  fdisc.f 
  setprob.f 
  user_4.6/qinit.f90 
  user_4.6/setaux.f90
  user_exceeds_th.f90
  ${all}/claw3_setaux_manifold.f90
  ${all}/rpn3_vc_advection.f90 
  ${all}/rpt3_vc_advection.f90 
  ${all}/rptt3_vc_advection.f90 
)

add_executable(filament3dx 
  filament.cpp 
  filament_options.c 
  filament_user.cpp 
  filament_map_extrude.c 
  ${all}/advection_patch_setup_manifold.c 
  ${all2d}/fclaw2d_map_cart.c 
  ${all2d}/fclaw2d_map_fivepatch.c 
  ${all2d}/fclaw2d_map_bilinear.c 
  $<TARGET_OBJECTS:filament3dx_f>
)
set_target_properties(filament3dx PROPERTIES OUTPUT_NAME filament)

target_include_directories(filament3dx PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(filament3dx PRIVATE 
 FORESTCLAW::CLAWPACK3_46
)

add_test(NAME clawpack_advection_3d_filament COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_3d_filament PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")