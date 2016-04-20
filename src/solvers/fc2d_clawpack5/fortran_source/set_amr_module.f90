SUBROUTINE set_amr_module(mwaves_in,mcapa_in,mthlim_in,method_in)
  USE amr_module, ONLY: mwaves, mcapa, method, mthlim, use_fwaves
  IMPLICIT NONE

  INTEGER, INTENT(in) :: mwaves_in, mcapa_in, method_in(7)
  INTEGER, INTENT(in) :: mthlim_in(mwaves_in)

  mwaves = mwaves_in
  mcapa = mcapa_in
  method = method_in
  mthlim = mthlim_in
  use_fwaves = .FALSE.
  write(*,*) "What is the mwaves", mwaves
  write(*,*) "What is the mcapa", mcapa
  write(*,*) "What is the mthlim_in", mthlim_in
  write(*,*) "What is the method_in", method_in
END SUBROUTINE set_amr_module
