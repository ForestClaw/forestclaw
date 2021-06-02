# -- fortran library
add_library(cudaclaw_f OBJECT
	fortran_source/cudaclaw_inlinelimiter.f
	fortran_source/cudaclaw_set_capacity.f
	fortran_source/cudaclaw_block.f
	fortran_source/cudaclaw_bc2_default.f
	fortran_source/cudaclaw_flux2.f
	fortran_source/cudaclaw_flux2fw.f
	fortran_source/cudaclaw_step2.f
	fortran_source/cudaclaw_step2_wrap.f
)

set_target_properties(cudaclaw_f PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

# -- add library
add_library(cudaclaw STATIC)
set_target_properties(cudaclaw PROPERTIES EXPORT_NAME CUDACLAW)

# device code needs to relcoateable since it is in a shared library
target_compile_options(cudaclaw PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
                       --relocatable-device-code=true
                       >)
# needed because we have device code in a library
set_target_properties(cudaclaw PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_sources(cudaclaw PRIVATE 
  fc2d_cudaclaw.cpp
	fc2d_cudaclaw_options.c
	cuda_source/cudaclaw_step2.cu
	cuda_source/cudaclaw_initialize.cu
	cuda_source/cudaclaw_limiters.cu
	cuda_source/cudaclaw_allocate.cu
	cuda_source/cudaclaw_flux2.cu
	fc2d_cuda_profiler.cu
	cuda_source/cudaclaw_store_patches.cu
  $<TARGET_OBJECTS:cudaclaw_f>
)

target_link_libraries(cudaclaw PUBLIC clawpatch CUDA::nvToolsExt)
target_include_directories(cudaclaw
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

# -- install
install(FILES
	fc2d_cudaclaw.h
	fc2d_cudaclaw_cuda.h
	fc2d_cudaclaw_options.h
	cudaclaw_user_fort.h
	fc2d_cudaclaw_fort.h
  DESTINATION include
)

install(TARGETS cudaclaw
  EXPORT ${PROJECT_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

# imported target, for use from FetchContent
add_library(FORESTCLAW::CUDACLAW INTERFACE IMPORTED GLOBAL)
target_link_libraries(FORESTCLAW::CUDACLAW INTERFACE cudaclaw)