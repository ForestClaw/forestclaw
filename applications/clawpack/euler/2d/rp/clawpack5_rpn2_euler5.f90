SUBROUTINE clawpack5_rpn2_euler5(ixy,maxm,meqn,mwaves,maux, &
   mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
!
! Roe-solver for the Euler equations with a tracer variable and separate shear
! and entropy waves.

! waves: 4
! equations: 5

! Conserved quantities:
!       1 density
!       2 x_momentum
!       3 y_momentum
!       4 energy
!       5 tracer

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.
! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
!     f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.  With the Roe solver we have
!      amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated into the
! flux differences.

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                   and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr
!
! This routine has been made thread safe by removing the common block storage
! of the Roe-averages.


IMPLICIT NONE

  ! Input
INTEGER, INTENT(in) :: ixy, maxm, meqn, mwaves, mbc, mx, maux
REAL(kind=8), INTENT(in) :: ql(meqn, 1-mbc:maxm+mbc)
REAL(kind=8), INTENT(in) :: qr(meqn, 1-mbc:maxm+mbc)
REAL(kind=8), INTENT(in) :: auxl(maux, 1-mbc:maxm+mbc)
REAL(kind=8), INTENT(in) :: auxr(maux, 1-mbc:maxm+mbc)

  ! Output
REAL(kind=8), INTENT(in out) :: wave(meqn, mwaves, 1-mbc:maxm+mbc)
REAL(kind=8), INTENT(in out) :: s(mwaves, 1-mbc:maxm+mbc)
REAL(kind=8), INTENT(in out) :: apdq(meqn, 1-mbc:maxm+mbc)
REAL(kind=8), INTENT(in out) :: amdq(meqn, 1-mbc:maxm+mbc)

  ! Local storage
INTEGER :: i, m, mw, mu, mv
REAL(kind=8) :: rho_sqrtl, rho_sqrtr, rho_sq2, u, v, u2v2, pl, pr, enth, c
REAL(kind=8) :: g1c2, euv
REAL(kind=8) :: delta(4), a(4)

REAL(kind=8) :: rhoim1, pim1, cim1, s0, rho1, rhou1, rhov1, en1, p1, c1, s1
REAL(kind=8) :: sfract, rhoi, pi, ci, s3, rho2, rhou2, rhov2, en2, p2, c2
REAL(kind=8) :: s2, df

  ! Use entropy fix for transonic rarefactions
LOGICAL, PARAMETER :: efix = .TRUE.

  ! Common block storage
  ! Ideal gas constant
REAL(kind=8) :: gamma, gamma1
COMMON /cparam/  gamma

gamma1 = gamma - 1.d0

  ! Set mu to point to  the component of the system that corresponds to
  ! momentum in the direction of this slice, mv to the orthogonal momentum:
IF (ixy == 1) THEN
   mu = 2
   mv = 3
ELSE
   mu = 3
   mv = 2
ENDIF
  ! Note that notation for u and v reflects assumption that the Riemann
  ! problems  are in the x-direction with u in the normal direciton and v in
  ! the orthogonal direcion, but with the above definitions of mu and mv the
  ! routine also works  with ixy=2 and returns, for example, f0 as the Godunov
  ! flux g0 for the Riemann problems u_t + g(u)_y = 0 in the y-direction.

  ! Initialize waves and speeds
wave = 0.d0
s = 0.d0

  ! Primary loop over grid cell interfaces
DO i = 2-mbc, mx+mbc
   !! Compute Roe-averaged quantities
   rho_sqrtl = SQRT(qr(1,i-1))
   rho_sqrtr = SQRT(ql(1,i))
   rho_sq2 = rho_sqrtl + rho_sqrtr

   u = (qr(mu,i-1) / rho_sqrtl + ql(mu,i) / rho_sqrtr) / rho_sq2
   v = (qr(mv,i-1) / rho_sqrtl + ql(mv,i) / rho_sqrtr) / rho_sq2
   u2v2 = u**2 + v**2

   pl = gamma1 * (qr(4,i-1) - 0.5d0 * (qr(2,i-1)**2 + qr(3,i-1)**2) / qr(1,i-1))
   pr = gamma1 * (ql(4,i)   - 0.5d0 * (ql(2,i)**2   + ql(3,i)**2) / ql(1,i))
   enth = (((qr(4,i-1) + pl) / rho_sqrtl + (ql(4,i) + pr) / rho_sqrtr)) / rho_sq2

   c = SQRT(gamma1 * (enth - 0.5d0 * u2v2))
   g1c2 = gamma1 / (gamma1 * (enth - 0.5d0 * u2v2))
   euv = enth - u2v2

      ! Now split the jump in q at each interface into waves and find a1 thru
      ! a4, the coefficients of the 4 eigenvectors:
   delta(1) = ql(1,i) - qr(1,i-1)
   delta(2) = ql(mu,i) - qr(mu,i-1)
   delta(3) = ql(mv,i) - qr(mv,i-1)
   delta(4) = ql(4,i) - qr(4,i-1)

   a(3) = g1c2 * (euv * delta(1) + u * delta(2) + v * delta(3) - delta(4))
   a(2) = delta(3) - v * delta(1)
   a(4) = (delta(2) + (c-u) * delta(1) - c * a(3)) / (2.d0 * c)
   a(1) = delta(1) - a(3) - a(4)

      ! Compute the waves
      ! Acoustic, left-going
   wave( 1,1,i) = a(1)
   wave(mu,1,i) = a(1) * (u - c)
   wave(mv,1,i) = a(1) * v
   wave( 4,1,i) = a(1) * (enth - u * c)
   wave( 5,1,i) = 0.d0

   s(1,i) = u - c

      ! Shear wave
   wave( 1,2,i) = 0.d0
   wave(mu,2,i) = 0.d0
   wave(mv,2,i) = a(2)
   wave( 4,2,i) = a(2) * v
   wave( 5,2,i) = 0.d0

   s(2,i) = u

      ! Entropy
   wave( 1,3,i) = a(3)
   wave(mu,3,i) = a(3) * u
   wave(mv,3,i) = a(3) * v
   wave( 4,3,i) = a(3) * 0.5d0 * u2v2
   wave( 5,3,i) = 0.d0

   s(3,i) = u

      ! Acoustic, right-going
   wave( 1,4,i) = a(4)
   wave(mu,4,i) = a(4) * (u + c)
   wave(mv,4,i) = a(4) * v
   wave( 4,4,i) = a(4) * (enth + u * c)
   wave( 5,4,i) = 0.d0

   s(4,i) = u + c

      ! Tracer wave - only carried in 5th component of q
   wave(5,5,i) = ql(5,i) - qr(5,i-1)

   s(5,i) = u

END DO

IF (.NOT.efix) THEN
      ! Compute flux differences amdq and apdq.
      ! No entropy fix
      !
      ! amdq = SUM s*wave   over left-going waves
      ! apdq = SUM s*wave   over right-going waves
   amdq = 0.d0
   apdq = 0.d0
   DO i=2-mbc, mx+mbc
      DO mw=1,mwaves
         IF (s(mw,i) < 0.d0) THEN
            amdq(:,i) = amdq(:,i) + s(mw,i) * wave(:,mw,i)
         ELSE
            apdq(:,i) = apdq(:,i) + s(mw,i) * wave(:,mw,i)
         ENDIF
      ENDDO
   ENDDO

ELSE
      ! With entropy fix

      ! compute flux differences amdq and apdq.
      ! First compute amdq as sum of s*wave for left going waves.
      ! Incorporate entropy fix by adding a modified fraction of wave
      ! if s should change sign.
   amdq = 0.d0
   apdq = 0.d0
   DO i = 2-mbc, mx+mbc
          ! Check 1-wave
      rhoim1 = qr(1,i-1)
      pim1 = gamma1*(qr(4,i-1) - 0.5d0*(qr(mu,i-1)**2 + qr(mv,i-1)**2) / rhoim1)
      cim1 = dsqrt(gamma*pim1/rhoim1)
      s0 = qr(mu,i-1)/rhoim1 - cim1     ! u-c in left state (cell i-1)

          ! Check for fully supersonic case:
      IF (s0 >= 0.d0 .AND. s(1,i) > 0.d0)  THEN
         !! Everything is right-going
         CYCLE
      ENDIF

      rho1 = qr(1,i-1) + wave(1,1,i)
      rhou1 = qr(mu,i-1) + wave(mu,1,i)
      rhov1 = qr(mv,i-1) + wave(mv,1,i)
      en1 = qr(4,i-1) + wave(4,1,i)
      p1 = gamma1 * (en1 - 0.5d0 * (rhou1**2 + rhov1**2) / rho1)
      c1 = SQRT(gamma * p1 / rho1)
      s1 = rhou1 / rho1 - c1  ! u-c to right of 1-wave
      IF (s0 < 0.d0 .AND. s1 > 0.d0) THEN
         !! Transonic rarefaction in the 1-wave
         sfract = s0 * (s1-s(1,i)) / (s1-s0)
      ELSE IF (s(1,i) < 0.d0) THEN
         !! 1-wave is leftgoing
         sfract = s(1,i)
      ELSE
         !! 1-wave is rightgoing
         sfract = 0.d0   ! This shouldn't happen since s0 < 0
      ENDIF
      DO m=1,meqn
         amdq(m,i) = sfract * wave(m,1,i)
      END DO

      !! Check contact discontinuity:
      IF (s(2,i) >= 0.d0) THEN
         !!# 2- 3- and 5-waves are rightgoing
         CYCLE
      ENDIF

      amdq(:,i) = amdq(:,i) + s(2,i) * wave(:,2,i)
      amdq(:,i) = amdq(:,i) + s(3,i) * wave(:,3,i)
      amdq(:,i) = amdq(:,i) + s(5,i) * wave(:,5,i)

      !! Check 4-wave:
      rhoi = ql(1,i)
      pi = gamma1*(ql(4,i) - 0.5d0*(ql(mu,i)**2 + ql(mv,i)**2) / rhoi)
      ci = SQRT(gamma*pi/rhoi)
      s3 = ql(mu,i)/rhoi + ci     ! u+c in right state  (cell i)

      rho2 = ql(1,i) - wave(1,4,i)
      rhou2 = ql(mu,i) - wave(mu,4,i)
      rhov2 = ql(mv,i) - wave(mv,4,i)
      en2 = ql(4,i) - wave(4,4,i)
      p2 = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2)/rho2)
      c2 = SQRT(gamma*p2/rho2)
      s2 = rhou2/rho2 + c2   ! u+c to left of 4-wave
      IF (s2 < 0.d0 .AND. s3 > 0.d0) THEN
         !! Transonic rarefaction in the 4-wave
         sfract = s2 * (s3-s(4,i)) / (s3-s2)
      ELSE IF (s(4,i) < 0.d0) THEN
         !! 4-wave is leftgoing
         sfract = s(4,i)
      ELSE
         !! 4-wave is rightgoing
         CYCLE
      ENDIF

      amdq(:,i) = amdq(:,i) + sfract * wave(:,4,i)
   ENDDO

   !! Compute the remaining right-going flux differences:
   !! df = SUM s*wave   is the total flux difference and apdq = df - amdq
   DO m=1,meqn
      DO i = 2-mbc, mx+mbc
         df = 0.d0
         DO mw=1,mwaves
            df = df + s(mw,i) * wave(m,mw,i)
         ENDDO
         apdq(m,i) = df - amdq(m,i)
      ENDDO
   ENDDO

ENDIF
END SUBROUTINE clawpack5_rpn2_euler5
