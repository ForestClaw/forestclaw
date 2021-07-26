subroutine cudaclaw_qinit(maxmx,maxmy,meqn,mbc,mx,my, & 
    xlower,ylower,dx,dy,q,maux,aux)

    !! # Set initial conditions for q.
    !! # Acoustics with smooth radially symmetric profile to test accuracy

    implicit none
    integer :: maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision :: xlower, ylower, dx, dy
    double precision :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
    double precision :: aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

    double precision :: pi, pi2
    common /compi/ pi, pi2

    integer i,j
    double precision :: xc, yc, rc, pressure, width

    width = 0.2d0

    do i = 1-mbc,mx+mbc
        xc = xlower + (i - 0.5d0)*dx
        do j = 1-mbc,my+mbc
            yc = ylower + (j - 0.5d0)*dy
            rc = dsqrt(xc**2 + yc**2)

            if (abs(rc-0.5d0) .le. width) then
                pressure = 1.d0 + cos(pi*(rc - 0.5d0)/width)
            else
                pressure = 0.d0
            endif
            q(i,j,1) = pressure
            q(i,j,2) = 0.d0
            q(i,j,3) = 0.d0
        end do
    end do
    return
end subroutine cudaclaw_qinit
