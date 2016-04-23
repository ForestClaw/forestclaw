SUBROUTINE geoclaw_set_amr_module(mwaves_in,mcapa_in,mthlim_in,method_in)
  USE amr_module, ONLY: mwaves, mcapa, method, mthlim, use_fwaves
  IMPLICIT NONE

  INTEGER, INTENT(in) :: mwaves_in, mcapa_in, method_in(7)
  INTEGER, INTENT(in) :: mthlim_in(mwaves_in)

  mwaves = mwaves_in
  mcapa = mcapa_in
  method = method_in
  mthlim = mthlim_in
  use_fwaves = .FALSE.
END SUBROUTINE geoclaw_set_amr_module
