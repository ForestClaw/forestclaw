SUBROUTINE clawpack5_rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr, &
     aux1,aux2,aux3,asdq,bmasdq,bpasdq)
  !!   =====================================================

  IMPLICIT DOUBLE PRECISION (a-h,o-z)

  !! # Riemann solver in the transverse direction for 2D Burgers' equation

  !! # Split asdq into eigenvectors of Roe matrix B.
  !! # For the scalar equation, this simply amounts to computing the
  !! # transverse wave speed from the opposite Riemann problem.

  DIMENSION    ql(meqn, 1-mbc:maxm+mbc)
  DIMENSION    qr(meqn, 1-mbc:maxm+mbc)
  DIMENSION   asdq(meqn, 1-mbc:maxm+mbc)
  DIMENSION bmasdq(meqn, 1-mbc:maxm+mbc)
  DIMENSION bpasdq(meqn, 1-mbc:maxm+mbc)

  !! # x- and y- Riemann problems are identical, so it doesn't matter if
  !! # ixy=1 or 2.

  DO i = 2-mbc, mx+mbc
     sb = 0.5d0*(qr(1,i-1) + ql(1,i))
     bmasdq(1,i) = dmin1(sb, 0.d0) * asdq(1,i)
     bpasdq(1,i) = dmax1(sb, 0.d0) * asdq(1,i)
  ENDDO

  RETURN
END SUBROUTINE clawpack5_rpt2
