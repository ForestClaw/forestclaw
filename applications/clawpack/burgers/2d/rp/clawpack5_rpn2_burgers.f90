!   ===============================================================
SUBROUTINE clawpack5_rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr, &
     auxl,auxr,wave,s,amdq,apdq)
!   ===============================================================

! Riemann solver for Burgers' equation in 2d:
!  u_t + 0.5*(u^2)_x + 0.5*(u^2)_y = 0

! waves: 1
! equations: 1

! Conserved quantities:
!       1 q

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q
!
! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr


        implicit double precision (a-h,o-z)

        dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
        dimension    s(mwaves, 1-mbc:maxm+mbc)
        dimension   ql(meqn, 1-mbc:maxm+mbc)
        dimension   qr(meqn, 1-mbc:maxm+mbc)
        dimension  apdq(meqn, 1-mbc:maxm+mbc)
        dimension  amdq(meqn, 1-mbc:maxm+mbc)
        logical efix

!     # x- and y- Riemann problems are identical, so it doesn't matter if
!     # ixy=1 or 2.

        efix = .true.

        DO  i = 2-mbc, mx+mbc
!        # wave is jump in q, speed comes from R-H condition:
            wave(1,1,i) = ql(1,i) - qr(1,i-1)
            s(1,i) = 0.5d0*(qr(1,i-1) + ql(1,i))

!        # compute left-going and right-going flux differences:

            amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
            apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)

            IF (efix) THEN
!               # entropy fix for transonic rarefactions:
                if (qr(1,i-1).lt.0.d0 .and. ql(1,i).gt.0.d0) then
                    amdq(1,i) = - 0.5d0 * qr(1,i-1)**2
                    apdq(1,i) =   0.5d0 * ql(1,i)**2
                end if
             END IF
          ENDDO

    return
  END SUBROUTINE clawpack5_rpn2
