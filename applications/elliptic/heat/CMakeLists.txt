add_library(heat_f OBJECT 
  hsmooth.f90 
  fortran/heat_setprob.f90 
  fortran/heat_fort_rhs.f90 
  fortran/heat_init.f90 
  fortran/heat_compute_error.f 
  fortran/heat_qexact.f90 
  fortran/heat_bc2.f 
  fortran/heat_apply_bc.f90 
  fortran/heat_eval_bc.f90 
  fortran/heat_tag4refinement.f 
  fortran/heat_tag4coarsening.f 
  fortran/heat_output.f90
)

add_library(heat_c OBJECT 
  heat_run.c 
  heat_options.c 
  heat_diagnostics.c 
)
target_link_libraries(heat_c PRIVATE 
  FORESTCLAW::FC2D_THUNDEREGG
)
add_executable(heat 
  heat.cpp 
  heat_user.cpp 
  $<TARGET_OBJECTS:heat_f>
  $<TARGET_OBJECTS:heat_c>
)

target_include_directories(heat PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(heat PRIVATE 
  FORESTCLAW::FC2D_THUNDEREGG ZLIB::ZLIB
)

add_test(NAME elliptic_heat 
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(elliptic_heat
         PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
