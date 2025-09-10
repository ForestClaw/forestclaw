SUBROUTINE fclaw3d_clawpatch46_fort_gauges_update(num,mx,my,mz,mbc,meqn,&
    xlower,ylower, zlower, dx,dy,dz,q,maux,aux,xc,yc,zc,qvar,avar)

    implicit none

    integer :: num, mx, my, mz,mbc, meqn, maux
    double precision :: xlower, ylower, zlower,dx, dy, dz, xc, yc, zc
    double precision :: q(1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc,maux)
    double precision :: qvar(meqn), avar(maux)

    !! local variables:
    double precision :: xcent,ycent,zcent,xoff,yoff,zoff
    integer :: iindex,jindex,kindex, mq, m, k

    double precision :: qz(0:1), az(0:1)

    iindex =  int((xc-xlower)/dx) + 1
    jindex =  int((yc-ylower)/dy) + 1
    kindex =  int((zc-zlower)/dz) + 1

    xcent  = xlower + (iindex-.5d0)*dx
    ycent  = ylower + (jindex-.5d0)*dy
    zcent  = zlower + (kindex-.5d0)*dz
    xoff   = (xc-xcent)/dx
    yoff   = (yc-ycent)/dy
    zoff   = (zc-zcent)/dz

    !! Linear interpolation between four cells
    do mq=1,meqn
        do k = 0,1
            qz(k) = (1.d0 - xoff) * (1.d0 - yoff) * q(iindex,jindex,kindex+k,mq)  &
            + xoff*(1.d0 - yoff) * q(iindex+1,jindex,kindex+k,mq)  &
            + (1.d0 - xoff) * yoff * q(iindex,jindex+1,kindex+k,mq)  &
            + xoff * yoff * q(iindex+1,jindex+1,kindex+k,mq)
        enddo
        qvar(mq) = qz(0) + (qz(1) - qz(0))*zoff
    enddo

    !! Linear interpolation between four cells
    do m=1,maux
        do k = 0,1
            az(k) = (1.d0 - xoff) * (1.d0 - yoff) *  aux(iindex,jindex,kindex+k,m)  &
            + xoff*(1.d0 - yoff) * aux(iindex+1,jindex,kindex+k,m)  &
            + (1.d0 - xoff) * yoff * aux(iindex,jindex+1,kindex+k,m)  &
            + xoff * yoff * aux(iindex+1,jindex+1,kindex+1,m)
        end do
        avar(m) = az(0) + (az(1) - az(0))*zoff
    enddo


    !! Zero out tiny values to prevent later problems reading data
    do mq = 1,meqn
        if (abs(qvar(mq)) < 1d-99) qvar(mq) = 0.d0
    end do
    do m = 1,maux
        if (abs(avar(m)) < 1d-99) avar(m) = 0.d0
    end do

END SUBROUTINE fclaw3d_clawpatch46_fort_gauges_update


!! This is a dummy routine that should be updated by the user.  
SUBROUTINE fclaw3d_clawpatch46_fort_gauges_move_local(num,mx,my,mz,mbc,meqn,&
    xlower,ylower, zlower, dx,dy,dz,q,maux,aux,xc,yc,zc, t, dt, &
    xc_new,yc_new,zc_new)

    implicit none

    integer :: num, mx, my, mz,mbc, meqn, maux
    double precision :: xlower, ylower, zlower,dx, dy, dz, xc, yc, zc
    double precision :: xc_new, yc_new, zc_new, t, dt
    double precision :: q(1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc,maux)

    xc_new = xc
    yc_new = yc
    zc_new = zc

END SUBROUTINE fclaw3d_clawpatch46_fort_gauges_move_local


