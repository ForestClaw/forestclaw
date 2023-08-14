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
            # remove everything including / before adding to the list
            list(APPEND non_empty_lines "${trimmed_line}")
        endif()
    endforeach()

    # Set the out_list to the filtered lines
    set(${out_list} "${non_empty_lines}" PARENT_SCOPE)
endfunction()

function(add_regressions filename)
    read_file_lines(${filename} tests) 

    set(working_directory "${CMAKE_CURRENT_SOURCE_DIR}")
    # Loop through each test
    foreach(test IN LISTS tests)
        string(REGEX REPLACE " " ";" argv "${test}")
        list(GET argv 0 executable)
        if(executable STREQUAL "cd")
            list(GET argv 1 directory)

            cmake_path(APPEND working_directory "${directory}")
            cmake_path(SET working_directory NORMALIZE "${working_directory}")

            message(STATUS "${working_directory}")
        else()
            string(REGEX REPLACE "^.*applications/" "" relative_source_directory "${working_directory}")
            set(test_name "${relative_source_directory}: ${test}")

            string(REGEX REPLACE "^\.*/" "" executable "${executable}")
            cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR "${executable}" OUTPUT_VARIABLE executable)

            list(REMOVE_AT argv 0)
            list(INSERT argv 0 "${executable}")


            message(STATUS "${argv}")

            # Add the test
            add_test(
                NAME 
                    ${test_name} 
                COMMAND 
                    ${MPIEXEC_EXECUTABLE} 
                    ${MPIEXEC_NUMPROC_FLAG}
                    ${MPIEXEC_MAX_NUMPROCS}
                    ${argv} 
                WORKING_DIRECTORY 
                    "${working_directory}" 
                COMMAND_EXPAND_LISTS
            )
        endif()

    endforeach()

endfunction()
