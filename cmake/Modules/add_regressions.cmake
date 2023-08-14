function(read_file_lines filepath out_list)
    # Check if the file exists
    if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${filepath})
        message(FATAL_ERROR "File not found: ${filepath}")
    endif()

    # Read the file content
    file(READ ${CMAKE_CURRENT_SOURCE_DIR}/${filepath} content)

    # Split the content by newline
    string(REPLACE "\n" ";" lines "${content}")

    set(non_empty_lines "")

    # Loop through each line
    foreach(line ${lines})
        # Trim whitespace from the beginning and end of the line
        string(STRIP "${line}" trimmed_line)

        # Check if the line is not empty and doesn't start with #
        if(NOT "${trimmed_line}" STREQUAL "" AND NOT "${trimmed_line}" MATCHES "^#.*")
            # Strip leading ./ before adding to the list
            string(REGEX REPLACE "^\\./" "" trimmed_line "${trimmed_line}")
            list(APPEND non_empty_lines "${trimmed_line}")
        endif()
    endforeach()

    # Set the out_list to the filtered lines
    set(${out_list} "${non_empty_lines}" PARENT_SCOPE)
endfunction()

function(add_regressions filename)
    read_file_lines(${filename} tests) 
    # Loop through each test
    foreach(test IN LISTS tests)
        set(test_command "${CMAKE_CURRENT_BINARY_DIR}/${test}")
        # esperate test into list by whitespace
        string(REGEX REPLACE " " ";" test_command "${test_command}")

        # remove everyting up to and including applications/ from source dir
        string(REGEX REPLACE "^.*applications/" "" test_name "${CMAKE_CURRENT_SOURCE_DIR}")
        set(test_name "${test_name}: ./${test}")

        message(STATUS "test         ${test}")
        message(STATUS "test_command ${test_command}")
        message(STATUS "test_args    ${args}")
        message(STATUS "test_name    ${test_name}")

        # Add the test
        add_test(NAME ${test_name} COMMAND ${test_command} ${test_args} WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" COMMAND_EXPAND_LISTS)

    endforeach()

endfunction()
