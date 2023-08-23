add_library(poisson_f OBJECT 
  hsmooth.f90 
  fortran/poisson_setprob.f90 
  fortran/poisson_fort_rhs.f90 
  fortran/poisson_compute_error.f 
  fortran/poisson_qexact.f90 
  fortran/poisson_apply_bc.f90 
  fortran/poisson_eval_bc.f90 
  fortran/poisson_output.f90
)

add_executable(poisson 
  poisson.cpp 
  poisson_user.cpp 
  poisson_options.c 
  poisson_diagnostics.c 
  $<TARGET_OBJECTS:poisson_f>
)

target_include_directories(poisson PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(poisson PRIVATE 
  FORESTCLAW::FC2D_THUNDEREGG
)

add_regressions(regressions.txt)
