set(all ${PROJECT_SOURCE_DIR}/clawpack/advection/p8/all)
set(all2d ${PROJECT_SOURCE_DIR}/clawpack/advection/2d/all)

add_executable(swirlp8
  swirl.cpp
  swirl_options.c
  ${all2d}/fclaw2d_map_cart.c
  ${all2d}/fclaw2d_map_fivepatch.c
  ${all2d}/fclaw2d_map_bilinear.c
)

set_target_properties(swirlp8 PROPERTIES OUTPUT_NAME swirl)

target_include_directories(swirlp8 PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${all}
)

target_link_libraries(swirlp8 PRIVATE
 FORESTCLAW::CLAWPACK3_46
)










