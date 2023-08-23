add_library(heat_phasefield_f OBJECT 
  heat/hsmooth.f90 
  heat/fortran/heat_setprob.f90 
  heat/fortran/heat_fort_rhs.f90 
  heat/fortran/heat_init.f90 
  heat/fortran/heat_compute_error.f 
  heat/fortran/heat_qexact.f90 
  heat/fortran/heat_bc2.f 
  heat/fortran/heat_apply_bc.f90 
  heat/fortran/heat_eval_bc.f90 
  heat/fortran/heat_tag4refinement.f 
  heat/fortran/heat_tag4coarsening.f 
  heat/fortran/heat_output.f90
  phasefield/cellave.f 
  phasefield/fortran/phasefield_setprob.f90 
  phasefield/fortran/phasefield_bc2.f 
  phasefield/fortran/phasefield_fort_rhs.f90 
  phasefield/fortran/phasefield_qinit.f90 
  phasefield/fortran/phasefield_apply_bc.f90 
  phasefield/fortran/phasefield_eval_bc.f90 
  phasefield/fortran/phasefield_tag4refinement.f 
  phasefield/fortran/phasefield_tag4coarsening.f
)

add_executable(heat_phasefield
  heat_phasefield.cpp 
  heat/heat.cpp 
  heat/heat_user.cpp 
  heat/heat_run.c 
  heat/heat_options.c 
  heat/heat_diagnostics.c 
  phasefield/phasefield.cpp
  phasefield/phasefield_operator.cpp 
  phasefield/phasefield_patch_operator.cpp 
  phasefield/phasefield_user.cpp 
  phasefield/phasefield_run.c 
  phasefield/phasefield_options.c 
  $<TARGET_OBJECTS:heat_phasefield_f>
)

target_include_directories(heat_phasefield PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(heat_phasefield PRIVATE 
  FORESTCLAW::FC2D_THUNDEREGG ZLIB::ZLIB
)

add_regressions(regressions.txt)
