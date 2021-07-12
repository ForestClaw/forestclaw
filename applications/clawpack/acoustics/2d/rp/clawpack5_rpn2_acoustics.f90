! =====================================================
SUBROUTINE clawpack5_rpn2(ixy,maxm,meqn,mwaves,maux,mbc,&
     mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================
! Riemann solver for the acoustics equations in 2d.

! waves: 2
! equations: 3

! Conserved quantities:
!       1 pressure
!       2 x_velocity
!       3 y_velocity

! Note that although there are 3 eigenvectors, the second eigenvalue
! is always zero and so we only need to compute 2 waves.
!
! solve Riemann problems along one slice of data.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q


! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr


  IMPLICIT DOUBLE PRECISION (a-h,o-z)

  DIMENSION wave(meqn, mwaves, 1-mbc:maxm+mbc)
  DIMENSION    s(mwaves, 1-mbc:maxm+mbc)
  DIMENSION   ql(meqn, 1-mbc:maxm+mbc)
  DIMENSION   qr(meqn, 1-mbc:maxm+mbc)
  DIMENSION apdq(meqn, 1-mbc:maxm+mbc)
  DIMENSION amdq(meqn, 1-mbc:maxm+mbc)
  DIMENSION auxl(maux, 1-mbc:maxm+mbc)
  DIMENSION auxr(maux, 1-mbc:maxm+mbc)

!     local arrays
!     ------------
  DIMENSION delta(3)

!     # density, bulk modulus, and sound speed, and impedence of medium:
!     # (should be set in setprob.f)
  COMMON /cparam/ rho,bulk,cc,zz



!     # set mu to point to  the component of the system that corresponds
!     # to velocity in the direction of this slice, mv to the orthogonal
!     # velocity:

  IF (ixy == 1) THEN
     mu = 2
     mv = 3
  ELSE
     mu = 3
     mv = 2
  ENDIF

!     # note that notation for u and v reflects assumption that the
!     # Riemann problems are in the x-direction with u in the normal
!     # direciton and v in the orthogonal direcion, but with the above
!     # definitions of mu and mv the routine also works with ixy=2
!     # in which case waves come from the
!     # Riemann problems u_t + g(u)_y = 0 in the y-direction.


!     # split the jump in q at each interface into waves

!     # find a1 and a2, the coefficients of the 2 eigenvectors:
  DO i = 2-mbc, mx+mbc
     delta(1) = ql(1,i) - qr(1,i-1)
     delta(2) = ql(mu,i) - qr(mu,i-1)
     a1 = (-delta(1) + zz*delta(2)) / (2.d0*zz)
     a2 = (delta(1) + zz*delta(2)) / (2.d0*zz)

    !        # Compute the waves.

     wave(1,1,i) = -a1*zz
     wave(mu,1,i) = a1
     wave(mv,1,i) = 0.d0
     s(1,i) = -cc

     wave(1,2,i) = a2*zz
     wave(mu,2,i) = a2
     wave(mv,2,i) = 0.d0
     s(2,i) = cc

  END DO


!     # compute the leftgoing and rightgoing flux differences:
!     # Note s(i,1) < 0   and   s(i,2) > 0.

  FORALL (m=1:meqn,  i=2-mbc: mx+mbc)
     amdq(m,i) = s(1,i)*wave(m,1,i)
     apdq(m,i) = s(2,i)*wave(m,2,i)
  END FORALL

  RETURN
END SUBROUTINE clawpack5_rpn2
