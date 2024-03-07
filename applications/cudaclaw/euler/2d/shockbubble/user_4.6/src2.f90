subroutine clawpack46_src2(maxmx,maxmy,meqn,mbc,mx,my, & 
           xlower,ylower,dx,dy,q,maux,aux,t,dt)
    !! =========================================================
    implicit none

    integer :: maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision :: xlower, ylower, dx, dy, t, dt
    double precision :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
    double precision :: aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

    !! # source terms for cylindrical symmetry in 2d Euler equations
    !! # about y = 0, so radius = y
    !! # y value at cell center is stored in aux(i,j,1)

    double precision gamma, gamma1
    common /cparam/  gamma,gamma1

    double precision qstar(4)
    double precision dt2, press, rad, rho, u, v

    integer :: i,j, ndim

    !! # 2-stage Runge-Kutta method

    dt2    = dt/2.d0
    press  = 0.d0
    ndim   = 2

    do  i = 1,mx
        do  j = 1,my
            rad      = ylower + (j-0.5d0)*dy
            rho      = q(i,j,1)
            u        = q(i,j,2)/q(i,j,1)
            v        = q(i,j,3)/q(i,j,1)
            press    = gamma1*(q(i,j,4) - 0.5d0*rho*(u**2 + v**2))

            if (rad .eq. 0.d0) then 
                write(6,*) 'src2.f90 : rad = 0'
            endif

            if (rho .eq. 0.d0) then 
                write(6,*) 'src2.f90 : rho = 0'
                stop
            endif

            qstar(1) = q(i,j,1) - dt2*(ndim-1)/rad * q(i,j,3)
            qstar(2) = q(i,j,2) - dt2*(ndim-1)/rad * (rho*u*v)
            qstar(3) = q(i,j,3) - dt2*(ndim-1)/rad * (rho*v*v)
            qstar(4) = q(i,j,4) - dt2*(ndim-1)/rad * v*(q(i,j,4) + press)

            !! # second stage

            rho      = qstar(1)
            u        = qstar(2)/qstar(1)
            v        = qstar(3)/qstar(1)
            press    = gamma1*(qstar(4) - 0.5d0*rho*(u**2 + v**2))
            if (rho.eq.0.d0) write(6,*) 'rho = 0 in qstar'

            q(i,j,1) = q(i,j,1) - dt*(ndim-1)/rad * qstar(3)
            q(i,j,2) = q(i,j,2) - dt*(ndim-1)/rad * (rho*u*v)
            q(i,j,3) = q(i,j,3) - dt*(ndim-1)/rad * (rho*v*v)
            q(i,j,4) = q(i,j,4) - dt*(ndim-1)/rad * v*(qstar(4) + press)
        end do
    end do
    return
end
