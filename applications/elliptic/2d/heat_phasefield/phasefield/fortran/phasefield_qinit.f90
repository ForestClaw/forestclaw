subroutine phasefield_init(meqn,mbc,mx,my, & 
        xlower,ylower,dx,dy,q)
    implicit none

    integer meqn,mbc,mx,my
    double precision xlower, ylower,dx,dy

    double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)

    integer i,j
    double precision xlow, ylow, w

    do j = 1-mbc,my+mbc
        do i = 1-mbc,mx+mbc
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy

            call cellave(xlow,ylow,dx,dy,w)
            q(i,j,1) = w-1.d0
            q(i,j,2) = 1-w
        enddo
    enddo

    return
end   subroutine phasefield_init


double precision function fdisc(x,y)
    implicit none

    double precision x,y

    DOUBLE PRECISION S, alpha, m, xi, k, gamma
    COMMON /comm_parms/ S, alpha, m, xi, k, gamma

    DOUBLE PRECISION r0, x0, y0
    COMMON /comm_init/ r0, x0, y0

    double precision r, rp, theta

    r = sqrt((x-x0)**2 + (y-y0)**2)
    theta = atan2(y - y0,x - x0)
    rp = r0*(1.d0 + gamma*cos(k*theta))
    fdisc = r - rp

end function fdisc

