SUBROUTINE fclaw2d_clawpatch5_gauge_update(num,mx,my,mbc,meqn,xlower,ylower, &
    dx,dy,q,maux,aux,xc,yc,qvar,avar)

    implicit none

    integer :: num, mx, my, mbc, meqn, maux
    double precision :: xlower, ylower, dx, dy, xc, yc
    double precision :: q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision :: qvar(meqn), avar(maux)

    !! local variables:
    double precision :: xcent,ycent,xoff,yoff
    integer :: iindex,jindex,mq, m

    iindex =  int((xc-xlower)/dx) + 1
    jindex =  int((yc-ylower)/dy) + 1

    xcent  = xlower + (iindex-.5d0)*dx
    ycent  = ylower + (jindex-.5d0)*dy
    xoff   = (xc-xcent)/dx
    yoff   = (yc-ycent)/dy

    !! Linear interpolation between four cells
    do mq=1,meqn
        qvar(mq) = (1.d0 - xoff) * (1.d0 - yoff) * q(mq,iindex,jindex)  &
        + xoff*(1.d0 - yoff) * q(mq,iindex+1,jindex)  &
        + (1.d0 - xoff) * yoff * q(mq,iindex,jindex+1)  &
        + xoff * yoff * q(mq,iindex+1,jindex+1)
    enddo

    !! Linear interpolation between four cells
    do m=1,maux
        avar(m) = (1.d0 - xoff) * (1.d0 - yoff) *  aux(m,iindex,jindex)  &
        + xoff*(1.d0 - yoff) * aux(m,iindex+1,jindex)  &
        + (1.d0 - xoff) * yoff * aux(m,iindex,jindex+1)  &
        + xoff * yoff * aux(m,iindex+1,jindex+1)
    enddo


    !! Zero out tiny values to prevent later problems reading data
    do mq = 1,meqn
        if (abs(qvar(mq)) < 1d-99) qvar(mq) = 0.d0
    end do
    do m = 1,maux
        if (abs(avar(m)) < 1d-99) avar(m) = 0.d0
    end do

END SUBROUTINE magic2d_update_gauge
