SUBROUTINE geoclaw_set_modules(mwaves_in,mcapa_in,mthlim_in,method_in)
  USE amr_module, ONLY: mwaves, mcapa, method, mthlim, use_fwaves
  USE geoclaw_module, ONLY: set_geo
!!  USE qinit_module, ONLY: set_qinit

  IMPLICIT NONE

  INTEGER, INTENT(in) :: mwaves_in, mcapa_in, method_in(7)
  INTEGER, INTENT(in) :: mthlim_in(mwaves_in)

!! Set values in amr_module
  mwaves = mwaves_in
  mcapa = mcapa_in
  method = method_in
  mthlim = mthlim_in
  use_fwaves = .FALSE.

!! Various modules from Geoclaw
  !! CALL set_geo()                    !# sets basic parameters g and coord system
!!  CALL set_qinit()                  !# specifies file with dh if this used instead

END SUBROUTINE geoclaw_set_modules
