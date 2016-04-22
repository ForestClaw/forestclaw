! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================
! Riemann-solver for the advection equation
!    q_t  +  u*q_x + v*q_y = 0
! where u and v are a given velocity field.

! waves: 1
! equations: 1
! aux fields: 2

! Conserved quantities:
!       1 q

! Auxiliary variables:
!         1  x_velocity
!         2  y_velocity

! solve Riemann problems along one slice of data.
! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! On output, wave contains the waves, s the speeds,
! and amdq, apdq the left-going and right-going flux differences,
! respectively.  Note that in this advective form, the sum of
! amdq and apdq is not equal to a difference of fluxes except in the
! case of constant velocities.

! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr

    implicit real*8(a-h,o-z)

    dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
    dimension    s(mwaves, 1-mbc:maxm+mbc)
    dimension   ql(meqn, 1-mbc:maxm+mbc)
    dimension   qr(meqn, 1-mbc:maxm+mbc)
    dimension apdq(meqn,1-mbc:maxm+mbc)
    dimension amdq(meqn,1-mbc:maxm+mbc)
    dimension auxl(maux,1-mbc:maxm+mbc)
    dimension auxr(maux,1-mbc:maxm+mbc)


!     # Set wave, speed, and flux differences:
!     ------------------------------------------

    do 30 i = 2-mbc, mx+mbc
        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        s(1,i) = auxl(ixy,i)
    !        # The flux difference df = s*wave  all goes in the downwind direction:
        amdq(1,i) = dmin1(auxl(ixy,i), 0.d0) * wave(1,1,i)
        apdq(1,i) = dmax1(auxl(ixy,i), 0.d0) * wave(1,1,i)
    30 END DO

    return
    end subroutine rpn2
