
add_library(testutils test.cpp)

if(TARGET ThunderEgg::ThunderEgg)
    target_sources(testutils PRIVATE utils/DomainReader.cpp)
endif(TARGET ThunderEgg::ThunderEgg)

target_include_directories(testutils PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(testutils PUBLIC ${PROJECT_BINARY_DIR})


target_include_directories(testutils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(testutils PUBLIC forestclaw)

if(TARGET ThunderEgg::ThunderEgg)
    target_link_libraries(testutils PUBLIC ThunderEgg::ThunderEgg)
endif(TARGET ThunderEgg::ThunderEgg)

if(mpi)
  target_link_libraries(testutils PUBLIC MPI::MPI_C)
endif(mpi)


macro(register_unit_tests test_executable)
  include(CTest)
  include(doctest)
  # relative path to executable
  string(REPLACE "${CMAKE_SOURCE_DIR}/" "" RELATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

  doctest_discover_tests(${test_executable} TEST_PREFIX "${RELATIVE_PATH}/${test_executable} --test-case\"" TEST_SUFFIX "\"" PROPERTIES LABELS "unit_test")
endmacro()