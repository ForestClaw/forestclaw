SUBROUTINE clawpack5_set_amr_module(mwaves_in,mcapa_in,mthlim_in, &
    method_in, use_fwaves_in)
  USE clawpack5_amr_module, ONLY: mwaves, mcapa, method, mthlim, use_fwaves
  IMPLICIT NONE

  INTEGER, INTENT(in) :: mwaves_in, mcapa_in, method_in(7)
  INTEGER, INTENT(in) :: mthlim_in(mwaves_in), use_fwaves_in

  integer :: mw

  mwaves = mwaves_in
  mcapa = mcapa_in
  method = method_in
  if (.not. allocated(mthlim)) then
      allocate(mthlim(mwaves))
  endif
  do mw = 1,mwaves
      mthlim(mw) = mthlim_in(mw)
  end do  
  use_fwaves = use_fwaves_in .ne. 0
END SUBROUTINE clawpack5_set_amr_module
