
add_library(square_f OBJECT 
  setprob.f
  q0.f
  fdisc.f
  cart_basis.f
  identity_basis.f
  mapping_terms.f
  square_velocity.f
  qexact.f
  dopri5.f
  user_4.6/square_output.f
  user_4.6/compute_error.f
  user_4.6/qinit.f
  user_4.6/square_setaux.f
  user_4.6/rpn2cons_fw_manifold.f
  user_4.6/rpt2cons_manifold.f
  user_4.6/rpn2cons_update_manifold.f
  user_4.6/tag4coarsening.f
  user_4.6/tag4refinement.f
)

add_executable(square 
  square.cpp
  square_user.cpp
  square_options.c
  fclaw2d_map_identity.c
  fclaw2d_map_cart.c
  fclaw2d_map_fivepatch.c
  fclaw2d_map_bilinear.c
  $<TARGET_OBJECTS:square_f>
)

target_include_directories(square PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
  ${advection_all}
)

target_link_libraries(square PRIVATE 
 FORESTCLAW::CLAWPACK4.6
 FORESTCLAW::CLAWPACK5
)