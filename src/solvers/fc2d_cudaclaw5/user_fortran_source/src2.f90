SUBROUTINE cudaclaw5_src2(meqn,mbc,mx,my, &
     xlower,ylower,dx,dy,q,maux,aux,t,dt)

  !! Called to update q by solving source term equation
  !! $q_t = \psi(q)$ over time dt starting at time t.
  !!
  !! This default version does nothing.

  IMPLICIT NONE
  INTEGER, INTENT(in) :: mbc,mx,my,meqn,maux
  REAL(kind=8), INTENT(in) :: xlower,ylower,dx,dy,t,dt
  REAL(kind=8), INTENT(in) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
  REAL(kind=8), INTENT(inout) ::  q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

END SUBROUTINE cudaclaw5_src2
