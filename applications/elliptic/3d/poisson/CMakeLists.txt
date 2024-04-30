add_library(poisson3d_f OBJECT 
  hsmooth.f90 
  fortran/poisson_setprob.f90 
  fortran/poisson_fort_rhs.f90 
  fortran/poisson_compute_error.f90 
  fortran/poisson_qexact.f90 
  fortran/poisson_apply_bc.f90 
  fortran/poisson_eval_bc.f90 
  fortran/poisson_output.f90
)

add_executable(poisson_3d 
  poisson.cpp 
  poisson_user.cpp 
  poisson_operator.cpp
  poisson_options.c 
  poisson_diagnostics.c 
  $<TARGET_OBJECTS:poisson3d_f>
)

# set name of executable
set_target_properties(poisson_3d
  PROPERTIES OUTPUT_NAME poisson
)

target_include_directories(poisson_3d PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(poisson_3d PRIVATE 
  FORESTCLAW::FC3D_THUNDEREGG
)

add_test(NAME elliptic_3d_poisson 
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(elliptic_3d_poisson
         PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
