set(rp ${PROJECT_SOURCE_DIR}/clawpack/euler/3d/rp)
set(all ${PROJECT_SOURCE_DIR}/clawpack/euler/3d/all)
set(all2d ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(overpressure_f OBJECT 
  setprob_module.f90 
  fdisc.f90
  latlong_exceeds_th.f90
  user_4.6/qinit.f90
  ${rp}/gethexinfo.f90
  ${rp}/rpn3_euler.f90
  ${rp}/rpt3_euler.f90
  ${rp}/rptt3_euler.f90
  ${all}/euler_setaux_manifold.f90
  ${rp}/rpn3_euler_mapped.f90
  ${rp}/rpt3_euler_mapped.f90
  ${rp}/rptt3_euler_mapped.f90
)

add_executable(overpressure 
  overpressure.cpp 
  overpressure_options.c 
  overpressure_user.cpp 
  ${all2d}/fclaw2d_map_cart.c 
  ${all2d}/fclaw2d_map_latlong.c 
  ${all2d}/fclaw2d_map_cubedsphere.c 
  ${all2d}/fclaw2d_map_pillowsphere.c 
  overpressure_map_extrude.c

  $<TARGET_OBJECTS:overpressure_f>
)

target_include_directories(overpressure PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(overpressure PRIVATE 
 FORESTCLAW::CLAWPACK3_46
)

add_test(NAME clawpack_euler_3d_overpressure COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_euler_3d_overpressure  PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")