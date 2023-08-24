add_library(allencahn_f OBJECT 
  hsmooth.f90 
  fortran/allencahn_setprob.f90 
  fortran/allencahn_bc2.f 
  fortran/allencahn_fort_rhs.f90 
  fortran/allencahn_qinit.f 
  fortran/allencahn_qexact.f90 
  fortran/allencahn_apply_bc.f90 
  fortran/allencahn_eval_bc.f90 
  fortran/allencahn_tag4refinement.f 
  fortran/allencahn_tag4coarsening.f
)

add_executable(allencahn 
  allencahn.cpp 
  allencahn_run.c 
  allencahn_user.cpp 
  allencahn_options.c 
  $<TARGET_OBJECTS:allencahn_f>
)

target_include_directories(allencahn PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
)

target_link_libraries(allencahn PRIVATE 
  FORESTCLAW::FC2D_THUNDEREGG
)

add_test(NAME elliptic_allencahn 
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regressions.sh
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(elliptic_allencahn
         PROPERTIES ENVIRONMENT "${FCLAW_TEST_ENVIRONMENT}")
