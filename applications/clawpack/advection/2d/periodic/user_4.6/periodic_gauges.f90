SUBROUTINE periodic_gauges_move_local(gnum,mx,my,mbc, & 
    meqn,xlower,ylower, dx,dy,q,maux,aux,xc,yc,t, dt, &
    xc_new,yc_new)

    implicit none

    integer :: gnum, mx, my, mbc, meqn, maux
    double precision :: xlower, ylower, dx, dy, xc, yc
    double precision :: xc_new,yc_new, t, dt
    double precision :: q(1-mbc:mx+mbc, 1-mbc:my+mbc,meqn)
    double precision :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc,maux)

    !! From setprob.data
    double precision uvel, vvel
    common /comvelocity/ uvel, vvel

    double precision pi, pi2
    common /compi/ pi,pi2

    !! Don't do anything yet ... See routine above for tips on how to get 
    !! information for the gauge

    if (gnum .eq. 0) then
        xc_new = 0.5*cos(pi*t)
        yc_new = 0.5*sin(pi*t)
    else if (gnum .lt. 10) then
        xc_new = xc + dt*uvel
        yc_new = yc + dt*vvel
    else        
        !! Gauges with ID > 10 do not move
        xc_new = xc
        yc_new = yc
    endif 


END SUBROUTINE periodic_gauges_move_local