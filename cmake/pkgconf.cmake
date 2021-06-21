# --- generate pkg-config .pc
set(pc_libs_private)
foreach(t clawpatch clawpack4.6 clawpack5 geoclaw)
  if(TARGET ${t})
    string(APPEND " -l${t}")
  endif()
endforeach()

set(pc_req_private "p4est sc ompi ompi-c orte")

set(pc_filename forestclaw.pc)
configure_file(${CMAKE_CURRENT_LIST_DIR}/pkgconf.pc.in ${pc_filename} @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${pc_filename} DESTINATION lib/pkgconfig)
