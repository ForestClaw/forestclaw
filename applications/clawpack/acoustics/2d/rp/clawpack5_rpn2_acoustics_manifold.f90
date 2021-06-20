! =====================================================
SUBROUTINE clawpack5_rpn2_manifold(ixy,maxm,meqn,mwaves, &
     maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================
! Riemann solver for the acoustics equations in 2d
! on general quadrilateral grid, with variable coefficients

! waves: 2
! equations: 3

! Conserved quantities:
!       1 pressure
!       2 x_velocity
!       3 y_velocity

! Auxiliary quantities:
!       1 a_x
!       2 a_y   where (a_x,a_y) is unit normal to left face
!       3 length_ratio_left ratio of length of left face to dyc
!       4 b_x
!       5 b_y   where (b_x,b_y) is unit normal to bottom face
!       6 length_ratio_bottom   ratio of length of bottom face to dxc
!       7 cell_area   ratio of cell area to dxc*dyc
!         (approximately Jacobian of mapping function)
!       8 Z (impedance)
!       9 c (sound speed)
!
! Rotate velocity and then call standard Riemann solver.
! The resulting waves and flux differences are then rotated
! back to x-y.

! solve Riemann problems along one slice of data.
!
! On input, ql contains the state vector at the left edge of each cell
! qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
! or the y-direction if ixy=2.
! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
! f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.
! With the Roe solver we have
!     amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated
! into the flux differences.
!
! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                   and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr.

  IMPLICIT NONE
!
  INTEGER, INTENT(in) :: maxm, meqn, mwaves, mbc, mx, maux
  DOUBLE PRECISION, INTENT(out) :: wave(meqn, mwaves, 1-mbc:maxm+mbc)
  DOUBLE PRECISION, INTENT(out) :: s(mwaves,  1-mbc:maxm+mbc)
  DOUBLE PRECISION, INTENT(in)  ::   ql(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION, INTENT(in)  ::   qr(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION, INTENT(out) :: apdq(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION, INTENT(out) :: amdq(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION, INTENT(in)  :: auxl(maux, 1-mbc:maxm+mbc)
  DOUBLE PRECISION, INTENT(in)  :: auxr(maux, 1-mbc:maxm+mbc)
!
!   ------------
  DOUBLE PRECISION :: delta(3)
  DOUBLE PRECISION :: unorl, unorr, a1, a2
  DOUBLE PRECISION :: alpha, beta, zi, zim, ci, cim
  INTEGER :: ixy, inx, iny, ilenrat, i, m
!
! Rotate the velocity vector (q(2), q(3)) so that it is aligned with the face
! normal.  The normal vector for the face at the i'th Riemann problem
! is stored in the aux array
! in locations (1,2) if ixy=1 or (4,5) if ixy=2.  The ratio of the
! length of the cell side to the length of the computational cell
! is stored in aux(3) or aux(6), respectively.

  IF (ixy.EQ.1) THEN
     inx = 1
     iny = 2
     ilenrat = 3
  ELSE
     inx = 4
     iny = 5
     ilenrat = 6
  ENDIF

! Determine rotation matrix:
!               [ alpha  beta ]
!               [-beta  alpha ]
! Note that this reduces to the identity on a standard cartesian grid.

! Determine normal velocity components at this edge:
  DO i=2-mbc,mx+mbc
     alpha = auxl(inx,i)
     beta  = auxl(iny,i)
     unorl = alpha*ql(2,i) + beta*ql(3,i)
     unorr = alpha*qr(2,i-1) + beta*qr(3,i-1)

     delta(1) = ql(1,i) - qr(1,i-1)
     delta(2) = unorl - unorr

     zi  = auxl(8,i)
     zim = auxl(8,i-1)
     ci  = auxl(9,i)
     cim = auxl(9,i-1)

     a1 = (-delta(1) + zi*delta(2)) / (zim + zi)
     a2 =  (delta(1) + zim*delta(2)) / (zim + zi)

!       Scale the velocities by the length ratios on each side.

     wave(1,1,i) = -a1*zim
     wave(2,1,i) = a1 * alpha
     wave(3,1,i) = a1 * beta
     s(1,i) = -cim * auxl(ilenrat,i)

     wave(1,2,i) = a2*zi
     wave(2,2,i) = a2 * alpha
     wave(3,2,i) = a2 * beta
     s(2,i) = ci * auxl(ilenrat,i)
  ENDDO

!   Compute the left-going and right-going fluctuations.
!   Note that s(1,i) < 0   and   s(2,i) > 0.

  DO m=1,meqn
     DO i = 2-mbc, mx+mbc
        amdq(m,i) = s(1,i)*wave(m,1,i)
        apdq(m,i) = s(2,i)*wave(m,2,i)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_rpn2_manifold
