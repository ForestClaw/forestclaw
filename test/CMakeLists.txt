
add_library(testutils test.cpp)

if(TARGET ThunderEgg::ThunderEgg)
    target_sources(testutils PRIVATE utils/DomainReader.cpp)
endif(TARGET ThunderEgg::ThunderEgg)

target_include_directories(testutils PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(testutils PUBLIC ${PROJECT_BINARY_DIR})


target_link_libraries(testutils PUBLIC forestclaw)

if(TARGET ThunderEgg::ThunderEgg)
    target_link_libraries(testutils PUBLIC ThunderEgg::ThunderEgg)
endif(TARGET ThunderEgg::ThunderEgg)

if(mpi)
  target_link_libraries(testutils PUBLIC MPI::MPI_C)
endif(mpi)

function(register_unit_tests test_target)
  include(CTest)
  include(Catch)
  catch_discover_tests(${test_target} TEST_PREFIX "unit_test: " PROPERTIES LABELS "unit_test")
endfunction()