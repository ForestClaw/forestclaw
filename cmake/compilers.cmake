include(CheckCCompilerFlag)

# --- compiler options

check_c_compiler_flag(-Wall _has_wall)
if(_has_wall)
  add_compile_options(
    $<$<COMPILE_LANGUAGE:C>:-Wall>
    $<$<COMPILE_LANGUAGE:CXX>:-Wall>
    $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
    )
else()
  check_c_compiler_flag(/Wall _has_msvc_wall)
  if(_has_msvc_wall)
    add_compile_options(/Wall)
  endif()
endif()

# check system networking libraries
include(CheckIncludeFile)
check_include_file(arpa/inet.h HAVE_ARPA_INET_H)
check_include_file(netinet/in.h HAVE_NETINET_IN_H)
if(WIN32 AND NOT HAVE_ARPA_INET_H AND NOT HAVE_NETINET_IN_H)
  check_include_file(Winsock2.h HAVE_WINSOCK2_H)
  set(WINSOCK_LIBRARIES wsock32 ws2_32 Iphlpapi)
endif()


# --- auto-ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()
