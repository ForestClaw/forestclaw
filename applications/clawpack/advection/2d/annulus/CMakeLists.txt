set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_library(annulus_f OBJECT 
  psi.f 
  setprob.f
  user_4.6/qinit.f 
  user_5.0/qinit.f90 
  ${all}/user46_setaux_manifold.f
  ${all}/user5_setaux_manifold.f90
  ${all}/clawpack46_rpn2adv_manifold.f 
  ${all}/clawpack46_rpt2adv_manifold.f 
  ${all}/clawpack5_rpn2adv_manifold.f90 
  ${all}/clawpack5_rpt2adv_manifold.f90
)

add_executable(annulus 
  annulus_user.cpp 
  annulus.cpp 
  annulus_options.c
	fclaw2d_map_annulus.c
  $<TARGET_OBJECTS:annulus_f>
)

target_include_directories(annulus PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(annulus PRIVATE 
  clawpack4.6
  clawpack5
)