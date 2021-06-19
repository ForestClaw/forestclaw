!   =====================================================
SUBROUTINE clawpack5_rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc, &
     mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
!   =====================================================
  IMPLICIT DOUBLE PRECISION (a-h,o-z)

!     # Riemann solver in the transverse direction for the acoustics equations.

!     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.

  DIMENSION     ql(meqn, 1-mbc:maxm+mbc)
  DIMENSION     qr(meqn, 1-mbc:maxm+mbc)
  DIMENSION   asdq(meqn, 1-mbc:maxm+mbc)
  DIMENSION bmasdq(meqn, 1-mbc:maxm+mbc)
  DIMENSION bpasdq(meqn, 1-mbc:maxm+mbc)

!     # density, bulk modulus, and sound speed, and impedence of medium:
!     # (should be set in setprob.f)
  COMMON /cparam/ rho,bulk,cc,zz


  IF (ixy == 1) THEN
     mu = 2
     mv = 3
  ELSE
     mu = 3
     mv = 2
  ENDIF

  DO  i = 2-mbc, mx+mbc
     a1 = (-asdq(1,i) + zz*asdq(mv,i)) / (2.d0*zz)
     a2 = (asdq(1,i) + zz*asdq(mv,i)) / (2.d0*zz)

    !        # The down-going flux difference bmasdq is the product  -c * wave

     bmasdq(1,i) = cc * a1*zz
     bmasdq(mu,i) = 0.d0
     bmasdq(mv,i) = -cc * a1

    !        # The up-going flux difference bpasdq is the product  c * wave

     bpasdq(1,i) = cc * a2*zz
     bpasdq(mu,i) = 0.d0
     bpasdq(mv,i) = cc * a2
  ENDDO

  RETURN
END SUBROUTINE clawpack5_rpt2
