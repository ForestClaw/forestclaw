SUBROUTINE geoclaw_set_modules(mwaves_in,mcapa_in,mthlim_in,method_in, ax, bx, ay, by)
  USE amr_module, ONLY: mwaves, mcapa, method, mthlim, use_fwaves, xlower, ylower, xupper, yupper

  !! These are not referenced below, but they compile.
  !!  USE regions_module, ONLY: set_regions
  !!  USE gauges_module, ONLY: set_gauges

  USE geoclaw_module, ONLY: set_geo
  USE topo_module, ONLY: read_topo_settings, read_dtopo_settings
  USE qinit_module, ONLY: set_qinit
  USE storm_module, only: set_storm
  USE friction_module, only: setup_variable_friction

  IMPLICIT NONE

  INTEGER, INTENT(in) :: mwaves_in, mcapa_in, method_in(7)
  INTEGER, INTENT(in) :: mthlim_in(mwaves_in)

  REAL(KIND=8), INTENT(IN) :: ax, bx, ay, by

!! Set values in amr_module
  mwaves = mwaves_in
  mcapa = mcapa_in
  method = method_in
  mthlim = mthlim_in
  use_fwaves = .FALSE.
  xlower = ax
  xupper = bx
  ylower = ay
  yupper = by

!! Various modules from Geoclaw
  CALL set_geo()                    !# sets basic parameters g and coord system
  CALL read_dtopo_settings()        !# specifies file with dtopo from earthquake
  CALL read_topo_settings()         !# specifies topography (bathymetry) files
  CALL set_qinit()                  !# specifies file with dh if this used instead
  CALL set_storm()                  ! Set storm parameters
  CALL setup_variable_friction()    ! Set variable friction parameters

END SUBROUTINE geoclaw_set_modules
