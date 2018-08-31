SUBROUTINE cudaclaw5_qinit(meqn,mbc,mx,my, &
     xlower,ylower,dx,dy,q,maux,aux)

  !!     # Set initial conditions for the q array.
  !!     # This default version prints an error message since it should
  !!     # not be used directly.  Copy this to an application directory and
  !!     # loop over all grid cells to set values of q(1:meqn, 1:mx, 1:my).

  IMPLICIT NONE

  INTEGER, INTENT(in) :: meqn,mbc,mx,my,maux
  REAL(kind=8), INTENT(in) :: xlower,ylower,dx,dy
  REAL(kind=8), INTENT(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
  REAL(kind=8), INTENT(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

  WRITE(6,*) '*** Error -- you must set initial conditions'
  STOP

END SUBROUTINE cudaclaw5_qinit
