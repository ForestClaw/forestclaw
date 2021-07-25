
set(rp ${PROJECT_SOURCE_DIR}/clawpack/shallow/2d/rp)

add_library(tsunami_f OBJECT
	setprob.f
	sgn1d/sgn.f90
	sgn1d/tri_solve.f90
	user_4.6/qinit.f90
	user_4.6/setaux.f90
	user_4.6/tag4refinement.f90
	user_4.6/tag4coarsening.f90
	bathy.f90
	${rp}/rpn2_tsunami.f90
	${rp}/rpt2_tsunami.f90
	${rp}/geoclaw_riemann_utils.f90
)

add_executable(tsunami
	tsunami.cpp
	tsunami_options.c
	tsunami_user.cpp
  $<TARGET_OBJECTS:tsunami_cuda_f>
)

target_include_directories(tsunami PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${rp}
)

target_link_libraries(tsunami PRIVATE 
	FORESTCLAW::CLAWPACK4.6
	FORESTCLAW::CLAWPACK5
)