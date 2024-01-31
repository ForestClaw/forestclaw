subroutine clawpack46_rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr, & 
    auxl,auxr,wave,s,amdq,apdq)
! Riemann-solver for the advection equation
!    q_t  +  u*q_x + v*q_y + w*q_z = 0
! where u and v are a given velocity field.
!
! waves: 1
! equations: 1
! aux fields: 3
!
! Conserved quantities:
!       1 q
!
! Auxiliary variables
!       1 x_velocity
!       2 y_velocity
!       3 z_velocity
!
! The velocities are specified at the left/bottom/back cell face centers.
!
! Solve Riemann problems along one slice of data.
! This data is along a slice in the x-direction if ixyz=1
!                               the y-direction if ixyz=2.
!                               the z-direction if ixyz=3.
!
! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
!
! On output, wave contains the waves, s the speeds,
! and amdq, apdq the left-going and right-going flux differences,
! respectively.  Note that in this advective form, the sum of
! amdq and apdq is not equal to a difference of fluxes except in the
! case of constant velocities.
!
! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr

    implicit none

    integer :: ixyz, maxm, meqn, mwaves, mbc, mx, maux

    double precision :: wave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision ::    s(mwaves,1-mbc:maxm+mbc)
    double precision ::   ql(meqn,1-mbc:maxm+mbc)
    double precision ::   qr(meqn,1-mbc:maxm+mbc)
    double precision :: amdq(meqn,1-mbc:maxm+mbc)
    double precision :: apdq(meqn,1-mbc:maxm+mbc)
    double precision :: auxl(maux,1-mbc:maxm+mbc)
    double precision :: auxr(maux,1-mbc:maxm+mbc)

    integer manifold
    common /com_manifold/ manifold

    integer :: i, iface, mcapa

    ! Set wave, speed, and flux differences:
    iface = ixyz
    mcapa = manifold
    do i = 2-mbc, mx+mbc
        wave(1,1,i) = ql(1,i) - qr(1,i-1)
        s(1,i) = auxl(iface+mcapa,i)
        ! The flux difference df = s*wave all goes in the downwind direction:
        amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
        apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)
    end do

    return
end subroutine clawpack46_rpn3
