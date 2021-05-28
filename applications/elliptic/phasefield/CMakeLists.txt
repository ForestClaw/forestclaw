add_library(phasefield_f OBJECT 
  cellave.f 
  fortran/phasefield_setprob.f90 
  fortran/phasefield_bc2.f 
  fortran/phasefield_fort_rhs.f90 
  fortran/phasefield_qinit.f90 
  fortran/phasefield_apply_bc.f90 
  fortran/phasefield_eval_bc.f90 
  fortran/phasefield_tag4refinement.f 
  fortran/phasefield_tag4coarsening.f
)

add_executable(phasefield 
  phasefield.cpp
  phasefield_run.c 
  phasefield_operator.cpp 
  phasefield_patch_operator.cpp 
  phasefield_user.cpp 
  phasefield_options.c 
  $<TARGET_OBJECTS:phasefield_f>
)

target_include_directories(phasefield PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(phasefield PRIVATE 
  FORESTCLAW::FC2D_THUNDEREGG
)

if(TARGET testutils)
  add_executable(phasefield_TEST.MPI1
    phasefield_run.c 
    phasefield_operator.cpp 
    phasefield_patch_operator.cpp 
    phasefield_patch_operator_TEST.cpp 
    phasefield_user.cpp 
    phasefield_options.c 
    $<TARGET_OBJECTS:phasefield_f>
  )
  
  target_link_libraries(phasefield_TEST.MPI1 PRIVATE 
    FORESTCLAW::FC2D_THUNDEREGG
    testutils
  )
endif(TARGET testutils)