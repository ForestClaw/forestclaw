SUBROUTINE fc2d_geoclaw_set_modules(mwaves_in, mcapa_in, meqn_in, maux_in,  & 
                                    mthlim_in, method_in, ax, bx, ay, by)
    USE amr_module, ONLY: mwaves, mcapa, method, mthlim, use_fwaves, xlower, ylower, xupper, yupper

    USE regions_module, ONLY: set_regions
!!    USE gauges_module, ONLY: set_gauges
    USE geoclaw_module, ONLY: set_geo
    USE topo_module, ONLY: read_topo_settings, read_dtopo_settings
    USE qinit_module, ONLY: set_qinit
    USE fixedgrids_module, ONLY: set_fixed_grids
    USE refinement_module, ONLY: set_refinement
    USE storm_module, only: set_storm
    USE friction_module, only: setup_variable_friction

    IMPLICIT NONE

    INTEGER, INTENT(in) :: mwaves_in, mcapa_in, method_in(7)
    INTEGER, INTENT(in) :: meqn_in, maux_in
    INTEGER, INTENT(in) :: mthlim_in(mwaves_in)

    !! We don't yet allow the user to specify a different gauges file
    CHARACTER(len=20) :: fname = 'gauges.data'

    INTEGER :: meqn, maux, mw

    REAL(KIND=8), INTENT(IN) :: ax, bx, ay, by

    logical :: restart !! This is read from claw.data;  set to .false. here

    !! Set values in amr_module
    mwaves = mwaves_in
    mcapa = mcapa_in
    meqn = meqn_in
    maux = maux_in
    method = method_in
    use_fwaves = .FALSE.
    xlower = ax
    xupper = bx
    ylower = ay
    yupper = by

    ALLOCATE(mthlim(mwaves))
    DO mw = 1,mwaves
       mthlim(mw) = mthlim_in(mw)
    ENDDO


    restart = .FALSE.

    !! Various modules from Geoclaw
    CALL set_geo()                    !# sets basic parameters g and coord system
    CALL set_regions()
    ! CALL set_gauges(restart,meqn,fname)
    CALL set_refinement()             !# sets refinement control parameters
    CALL read_dtopo_settings()        !# specifies file with dtopo from earthquake
    CALL read_topo_settings()         !# specifies topography (bathymetry) files
    CALL set_qinit()                  !# specifies file with dh if this used instead
    CALL set_fixed_grids()            !# Fixed grid settings
    CALL set_storm()                  ! Set storm parameters
    CALL setup_variable_friction()    ! Set variable friction parameters

END SUBROUTINE fc2d_geoclaw_set_modules
