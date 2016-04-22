SUBROUTINE geoclaw_setup()

  USE geoclaw_module, ONLY: set_geo
  USE topo_module, ONLY: read_topo_settings, read_dtopo_settings
  USE qinit_module, ONLY: set_qinit

  IMPLICIT NONE

  CALL set_geo()                    !# sets basic parameters g and coord system
  CALL read_dtopo_settings()        !# specifies file with dtopo from earthquake
  CALL read_topo_settings()         !# specifies topography (bathymetry) files
  CALL set_qinit()                  !# specifies file with dh if this used instead

  END SUBROUTINE geoclaw_setup
