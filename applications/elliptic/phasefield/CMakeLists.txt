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

add_test(NAME elliptic_phasefield 
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(elliptic_phasefield
         PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
