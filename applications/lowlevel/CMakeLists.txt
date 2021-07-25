# Test the low-level logic for smooth refinement
add_executable(smooth smooth.c)
target_link_libraries(smooth PRIVATE FORESTCLAW::FORESTCLAW)
