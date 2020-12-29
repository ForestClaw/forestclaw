subroutine clawpack46_src2(maxmx,maxmy,meqn,mbc,mx,my, & 
    xlower,ylower,dx,dy,q,maux,aux,t,dt)
    implicit none

    integer maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision xlower, ylower, dx, dy, t, dt
    double precision   q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
    double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

    integer i,j
    double precision xc,yc, xm(2), tm, delta, r, dsrc, tf

    xm(1) = 0.21d0
    xm(2) = 0.5d0

    tm = 0.36d0
    tf = 0.5d0

    do i = 1-mbc,mx+mbc
        xc = xlower + (i-0.5d0)*dx
        do j = 1-mbc,my+mbc
            yc = ylower + (j-0.5d0)*dy
            r = sqrt((xc-xm(1))**2 + (yc-xm(2))**2)  
            dsrc = delta(r)*delta(tf-tm-t)
            q(i,j,1) = q(i,j,1) + dt*dsrc
            end do
    end do
      

!!
    return
end subroutine clawpack46_src2
