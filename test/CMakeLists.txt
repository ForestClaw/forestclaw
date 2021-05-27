
add_library(testutils test.cpp utils/DomainReader.cpp)
target_include_directories(testutils PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(testutils PUBLIC ${PROJECT_BINARY_DIR})
target_link_libraries(testutils PUBLIC forestclaw)
