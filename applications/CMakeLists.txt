# --------------------------------------------
# List of available applications and examples
#
# Note : Any files below which are uncommented
# need to be available at compile time, even 
# the configuration flag is not set.  
# --------------------------------------------


cmake_minimum_required(VERSION 3.13...3.20)

project(ForestClawApps LANGUAGES C CXX Fortran)

# -- check if this being included as a subdirectory
if(NOT TARGET FORESTCLAW::FORESTCLAW)

    find_package(MPI)

    # --- get installed ForestClaw package
    find_package(FORESTCLAW REQUIRED)

endif()

# --------------------------------------------
# Clawpack 4.6 and 5.0 examples
# --------------------------------------------

if(TARGET FORESTCLAW::CLAWPACK4.6 AND TARGET FORESTCLAW::CLAWPACK5)

    include(clawpack/clawpack.cmake)

endif()

# --------------------------------------------
# ThunderEgg solver
# --------------------------------------------

if(TARGET FORESTCLAW::FC2D_THUNDEREGG)

    include(elliptic/elliptic.cmake)

endif()




# ----------------------------------
# Projects under development
# ----------------------------------
if(TARGET FORESTCLAW::CLAWPACK4.6 AND TARGET FORESTCLAW::CLAWPACK5)

    # Development projects
    #include(development.cmake)

endif()

# ----------------------------------
# Projects for paper
# ----------------------------------
if(TARGET FORESTCLAW::CLAWPACK4.6)

    #include(paper.cmake)

endif()


# ----------------------------------
# GeoClaw library and examples
# ----------------------------------

if(TARGET FORESTCLAW::GEOCLAW)

    include(geoclaw/geoclaw.cmake)

endif()
