set(rp ${PROJECT_SOURCE_DIR}/clawpack/burgers/2d/rp)

add_library(burgers_f OBJECT 
	user_4.6/qinit.f
	user_4.6/setaux.f
	user_4.6/b4step2.f
	user_4.6/limiter.f
	user_4.6/rpn2bu.f
	user_4.6/rpt2bu.f
	user_4.6/flux2.f
	user_4.6/rpn2cons_update.f
	user_4.6/rpn2cons_update_order2.f
	${rp}/clawpack46_rpn2_burgers.f 
	${rp}/clawpack46_rpt2_burgers.f 
	user_4.6/tag4refinement.f
	user_4.6/tag4coarsening.f
	user_4.6/periodic_interpolate.f
)

add_executable(burgers 
	burgers.cpp
	burgers_user.cpp
	burgers_options.c
  $<TARGET_OBJECTS:burgers_f>
)

target_include_directories(burgers PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(burgers PRIVATE 
  FORESTCLAW::CLAWPACK4.6
  FORESTCLAW::CLAWPACK5
)