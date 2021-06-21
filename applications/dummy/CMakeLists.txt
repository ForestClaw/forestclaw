add_executable(dummy 
    fclaw2d_dummy.c
    dummy_blackbox.c
)

target_include_directories(dummy PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${all}
)

target_link_libraries(dummy PRIVATE 
 FORESTCLAW::CLAWPATCH
)