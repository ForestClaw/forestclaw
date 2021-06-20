! =====================================================
SUBROUTINE clawpack5_rpn2_acoustics_vc(ixy,maxm,meqn,mwaves,maux,mbc,mx,&
     ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! Riemann solver for the acoustics equations in 2d, with varying
! material properties rho and kappa

! waves: 2
! equations: 3
! aux fields: 2

! Conserved quantities:
!       1 pressure
!       2 x_momentum
!       3 y_momentum

! Auxiliary variables:
!         1  density
!         2  sound_speed

! Note that although there are 3 eigenvectors, the second eigenvalue
! is always zero and so we only need to compute 2 waves.
!
! solve Riemann problems along one slice of data.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! Here it is assumed that auxl=auxr gives the cell values.

! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr


  IMPLICIT none

  integer  :: ixy,maxm,meqn,mwaves,maux,mbc,mx
  DOUBLE PRECISION :: wave(meqn, mwaves, 1-mbc:maxm+mbc)
  DOUBLE PRECISION ::    s(mwaves, 1-mbc:maxm+mbc)
  DOUBLE PRECISION ::   ql(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION ::   qr(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION :: apdq(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION :: amdq(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION :: auxl(maux, 1-mbc:maxm+mbc)
  DOUBLE PRECISION :: auxr(maux, 1-mbc:maxm+mbc)

!     local arrays
!     ------------
  integer i, mu, mv, m
  DOUBLE PRECISION :: delta(3), a1, a2, zi, zim

!     # set mu to point to  the component of the system that corresponds
!     # to velocity in the direction of this slice, mv to the orthogonal
!     # velocity.


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


!     # split the jump in q at each interface into waves
!     # The jump is split into a leftgoing wave traveling at speed -c
!     # relative to the material properties to the left of the interface,
!     # and a rightgoing wave traveling at speed +c
!     # relative to the material properties to the right of the interface,

!     # find a1 and a2, the coefficients of the 2 eigenvectors:
  DO  i = 2-mbc, mx+mbc
     delta(1) = ql(1,i) - qr(1,i-1)
     delta(2) = ql(mu,i) - qr(mu,i-1)
    !        # impedances:
     zi = auxl(1,i)*auxl(2,i)
     zim = auxl(1,i-1)*auxl(2,i-1)

     a1 = (-delta(1) + zi*delta(2)) / (zim + zi)
     a2 =  (delta(1) + zim*delta(2)) / (zim + zi)


    !        # Compute the waves.

     wave(1,1,i) = -a1*zim
     wave(mu,1,i) = a1
     wave(mv,1,i) = 0.d0
     s(1,i) = -auxl(2,i-1)

     wave(1,2,i) = a2*zi
     wave(mu,2,i) = a2
     wave(mv,2,i) = 0.d0
     s(2,i) = auxl(2,i)

  END DO



!     # compute the leftgoing and rightgoing flux differences:
!     # Note s(1,i) < 0   and   s(2,i) > 0.

  DO  m=1,meqn
     DO  i = 2-mbc, mx+mbc
        amdq(m,i) = s(1,i)*wave(m,1,i)
        apdq(m,i) = s(2,i)*wave(m,2,i)
     END DO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_rpn2_acoustics_vc
