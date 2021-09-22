set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/3d/all)

add_library(filament_f OBJECT 
  psi.f 
  fdisc.f 
  setprob.f 
  user_4.6/qinit.f90 
  user_4.6/setaux.f90
  ${all}/rpn3_vc_advection.f90 
  ${all}/rpt3_vc_advection.f90 
  ${all}/rptt3_vc_advection.f90 
)

add_executable(filament 
  filament.cpp 
  filament_options.c 
  filament_user.cpp 
  $<TARGET_OBJECTS:filament_f>
)

target_include_directories(filament PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(filament PRIVATE 
 FORESTCLAW::CLAWPACK3_46
)

add_test(NAME clawpack_advection_3d_filament COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(clawpack_advection_3d_filament PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")