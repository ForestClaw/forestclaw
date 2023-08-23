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

function(add_regressions filename source_dir build_dir mpiexec mpiexec_np_flag mpiexec_max_np)
    read_file_lines(${filename} tests) 

    set(working_directory "${source_dir}")
    # Loop through each test
    foreach(test IN LISTS tests)
            # Initialize an empty list to hold the final split values
            set(argv "")

            # Initialize a variable to hold the current token
            set(current_token "")

            # Initialize a variable to indicate if we are inside quotes
            set(in_quotes FALSE)

            # Get length of the string
            string(LENGTH ${test} test_length)

            # Loop through each character in the string
            foreach(index RANGE 0 ${test_length})
                # Use string(SUBSTRING) to get a single character from the string
                string(SUBSTRING ${test} ${index} 1 char)
              if(char STREQUAL "\"" OR char STREQUAL "'")
                # Toggle the in-quotes state
                set(in_quotes NOT ${in_quotes})
              elseif((char STREQUAL " " OR char STREQUAL "\t") AND NOT ${in_quotes})
                # If the character is a space or a tab and we are not in quotes,
                # then this is a delimiter.
                # Add the current token to the list and clear it
                if (NOT "${current_token}" STREQUAL "")
                  list(APPEND argv "${current_token}")
                endif()
                set(current_token "")
              else()
                # Otherwise, append the character to the current token
                set(current_token "${current_token}${char}")
              endif()
            endforeach()

            # Append the last token to the list
            if (NOT "${current_token}" STREQUAL "")
              list(APPEND argv "${current_token}")
            endif()


        list(GET argv 0 executable)
        list(REMOVE_AT argv 0)
        string(REGEX REPLACE "^.*/" "" executable "${executable}")

        if(executable STREQUAL "cd")
            list(GET argv 1 directory)

            cmake_path(APPEND working_directory "${directory}")
            cmake_path(SET working_directory NORMALIZE "${working_directory}")

            message(STATUS "${working_directory}")
        else()
            string(REGEX REPLACE "^.*applications/" "" relative_source_directory "${working_directory}")
            set(test_name "${relative_source_directory}: ${test}")

            cmake_path(APPEND build_dir "${executable}" OUTPUT_VARIABLE executable)

              add_test(${test_name} ${mpiexec} ${mpiexec_np_flag} ${mpiexec_max_np} ${executable} ${argv})

            set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY ${working_directory})

        endif()

    endforeach()

endfunction()
