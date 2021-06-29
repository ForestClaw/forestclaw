# --------------------------------------------
# List of available applications and examples
#
# Note : Any files below which are uncommented
# need to be available at compile time, even 
# the configuration flag is not set.  
# --------------------------------------------


cmake_minimum_required(VERSION 3.13...3.20)

project(ForestClawApps LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

# -- check if this being included as a subdirectory
if(NOT TARGET FORESTCLAW::FORESTCLAW)

    set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
    find_package(MPI)

    # --- get installed ForestClaw package
    find_package(FORESTCLAW REQUIRED)

    if(TARGET FORESTCLAW::CUDACLAW)
        enable_language(CUDA)
    endif()

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
if(TARGET FORESTCLAW::CLAWPACK4.6 AND TARGET FORESTCLAW::CLAWPACK5)

    include(paper/paper.cmake)

endif()


# ----------------------------------
# GeoClaw library and examples
# ----------------------------------

if(TARGET FORESTCLAW::GEOCLAW)

    include(geoclaw/geoclaw.cmake)

endif()

## ---------------------------------------------
## Ash3d examples (requires submodule)
## -- Uncomment Makefiles below to compile Ash3d
## ---------------------------------------------
#
#if FCLAW_ENABLE_ASH3D
#
### Solver library
##include src/solvers/fc3d_ash3d/Makefile.am
#
### Applications
##include src/solvers/fc3d_ash3d/test_cases/ash3d.apps
#
#endif
#
## ----------------------------------------
## ACTIVEFLUX project (requires submodule)
## ----------------------------------------
#if FCLAW_ENABLE_ACTIVEFLUX
#
### Patches needed for Active Flux solver
##include src/patches/clawpatch/Makefile.am
#
### Solver library
##include src/solvers/fc2d_activeflux/Makefile.am
#
### Applications
##include src/solvers/fc2d_activeflux/examples/activeflux.apps
#
#endif


# ----------------------------------
# CudaClaw library and examples
# ----------------------------------
if(TARGET FORESTCLAW::CUDACLAW)

    include(cudaclaw.cmake)

endif()


# Include examples that use the low-level interface to p4est
# They should always compile

add_subdirectory(lowlevel)