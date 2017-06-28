SUBROUTINE clawpack5_rpt2_euler5(ixy,imp,maxm,meqn,mwaves, &
     maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)

  IMPLICIT NONE

    ! Input
  INTEGER, INTENT(in) :: ixy, imp, maxm, meqn, mwaves, mbc, mx, maux
  REAL(kind=8), INTENT(in) :: ql(meqn, 1-mbc:maxm+mbc)
  REAL(kind=8), INTENT(in) :: qr(meqn, 1-mbc:maxm+mbc)
  REAL(kind=8), INTENT(in) :: aux1(maux, 1-mbc:maxm+mbc)
  REAL(kind=8), INTENT(in) :: aux2(maux, 1-mbc:maxm+mbc)
  REAL(kind=8), INTENT(in) :: aux3(maux, 1-mbc:maxm+mbc)
  REAL(kind=8), INTENT(in) :: asdq(meqn, 1-mbc:maxm+mbc)

    ! Output
  REAL(kind=8), INTENT(in out) :: bmasdq(meqn, 1-mbc:maxm+mbc)
  REAL(kind=8), INTENT(in out) :: bpasdq(meqn, 1-mbc:maxm+mbc)

    ! Local storage
  INTEGER :: i, mw, mu, mv
  REAL(kind=8) :: a(4), waveb(5,4), sb(4), rho_sqrtl, rho_sqrtr, rho_sq2
  REAL(kind=8) :: u, v, u2v2, pl, pr, enth, c, g1c2, euv

    ! Common block - ideal gas constant
  REAL(kind=8) :: gamma, gamma1
  COMMON /cparam/  gamma

  gamma1 = gamma - 1.d0

    ! Initialize output
  bmasdq = 0.d0
  bpasdq = 0.d0

    ! Normal and transverse direction indices
  IF (ixy == 1) THEN
     mu = 2
     mv = 3
  ELSE
     mu = 3
     mv = 2
  ENDIF

    ! Primary loop over grid cell interfaces
  DO i = 2-mbc, mx+mbc
     !! Roe averaged variables
     rho_sqrtl = dsqrt(qr(1,i-1))
     rho_sqrtr = dsqrt(ql(1,i))
     rho_sq2 = rho_sqrtl + rho_sqrtr
     u = (qr(mu,i-1) / rho_sqrtl + ql(mu,i) / rho_sqrtr) / rho_sq2
     v = (qr(mv,i-1) / rho_sqrtl + ql(mv,i) / rho_sqrtr) / rho_sq2
     u2v2 = u**2 + v**2

     pl = gamma1*(qr(4,i-1) - 0.5d0*(qr(2,i-1)**2 + qr(3,i-1)**2)/qr(1,i-1))
     pr = gamma1*(ql(4,i)   - 0.5d0*(ql(2,i)**2   + ql(3,i)**2)  /ql(1,i))
     enth = (((qr(4,i-1)+pl)/rho_sqrtl + (ql(4,i)+pr)/rho_sqrtr)) / rho_sq2

     c = SQRT(gamma1*(enth - 0.5d0*u2v2))
     g1c2 = gamma1 / c**2
     euv = enth - u2v2

     a(3) = g1c2 * (euv*asdq(1,i) + u*asdq(mu,i) + v*asdq(mv,i) - asdq(4,i))
     a(2) = asdq(mu,i) - u*asdq(1,i)
     a(4) = (asdq(mv,i) + (c-v)*asdq(1,i) - c*a(3)) / (2.d0*c)
     a(1) = asdq(1,i) - a(3) - a(4)

     waveb(1,1) = a(1)
     waveb(mu,1) = a(1)*u
     waveb(mv,1) = a(1)*(v-c)
     waveb(4,1) = a(1)*(enth - v*c)
     waveb(5,1) = 0.d0
     sb(1) = v - c

     waveb(1,2) = a(3)
     waveb(mu,2) = a(3)*u + a(2)
     waveb(mv,2) = a(3)*v
     waveb(4,2) = a(3)*0.5d0*u2v2 + a(2)*u
     waveb(5,2) = 0.d0
     sb(2) = v

     waveb(1,3) = a(4)
     waveb(mu,3) = a(4)*u
     waveb(mv,3) = a(4)*(v+c)
     waveb(4,3) = a(4)*(enth+v*c)
     waveb(5,3) = 0.d0
     sb(3) = v + c

     waveb(1,4) = 0.d0
     waveb(mu,4) = 0.d0
     waveb(mv,4) = 0.d0
     waveb(4,4) = 0.d0
     waveb(5,4) = asdq(5,i)
     sb(4) = v

        ! Compute the flux differences bmasdq and bpasdq
     DO mw=1,4
        bmasdq(:,i) = bmasdq(:,i) + MIN(sb(mw), 0.d0) * waveb(:,mw)
        bpasdq(:,i) = bpasdq(:,i) + MAX(sb(mw), 0.d0) * waveb(:,mw)
     END DO
  END DO

END SUBROUTINE clawpack5_rpt2_euler5
