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

add_executable(heat 
  heat.cpp 
  heat_run.c 
  heat_user.cpp 
  heat_options.c 
  heat_diagnostics.c 
  $<TARGET_OBJECTS:heat_f>
)

target_include_directories(heat PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(heat PRIVATE 
  thunderegg
)