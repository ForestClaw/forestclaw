# -- fortran library
add_library(cudaclaw5_f OBJECT
	fortran_source/cudaclaw5_amr_module.f90
	fortran_source/cudaclaw5_set_amr_module.f90
	fortran_source/cudaclaw5_inlinelimiter.f
	fortran_source/cudaclaw5_limiter.f
	fortran_source/cudaclaw5_set_capacity.f
	fortran_source/cudaclaw5_block.f
	fortran_source/cudaclaw5_bc2_default.f
	fortran_source/cudaclaw5_flux2.f
	fortran_source/cudaclaw5_step2.f90
	fortran_source/cudaclaw5_step2_wrap.f
	fortran_source/cudaclaw5_update_q.f
)

set_target_properties(cudaclaw5_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(cudaclaw5 STATIC)
set_target_properties(cudaclaw5 PROPERTIES EXPORT_NAME CUDACLAW5)

# device code needs to relcoateable since it is in a shared library
target_compile_options(cudaclaw5 PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
                       --relocatable-device-code=true
                       >)
# needed because we have device code in a library
set_target_properties(cudaclaw5 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_sources(cudaclaw5 PRIVATE 
	fc2d_cudaclaw5.cpp
	fc2d_cudaclaw5_options.c
	cuda_source/cudaclaw5_update_q.cu
	cuda_source/cudaclaw5_step2.cu
	cuda_source/cudaclaw5_allocate.cu
	cuda_source/cudaclaw5_flux2.cu
	$<TARGET_OBJECTS:cudaclaw_f>
)

target_link_libraries(cudaclaw5 PUBLIC clawpatch CUDA::nvToolsExt)
target_include_directories(cudaclaw5
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

# -- install
install(FILES
	fc2d_cudaclaw5.h
	fc2d_cudaclaw5_options.h
	cudaclaw5_user_fort.h
	fc2d_cudaclaw5_fort.h
	cuda_source/cudaclaw5_update_q.h
	cuda_source/cudaclaw5_allocate.h
	DESTINATION include
)

install(TARGETS cudaclaw5
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

# imported target, for use from FetchContent
add_library(FORESTCLAW::CUDACLAW5 INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::CUDACLAW5 INTERFACE cudaclaw5)