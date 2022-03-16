SUBROUTINE fc2d_geoclaw_update_gauge (mx,my,mbc,meqn,xlower,ylower, &
    dx,dy,q,maux,aux,xc,yc,qvar,avar)

    use geoclaw_module, only: dry_tolerance

    implicit none

    integer :: mx, my, mbc, meqn, maux
    real(kind=8) :: xlower, ylower, dx, dy, xc, yc
    real(kind=8) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8) :: qvar(meqn), avar(maux)


    !! local variables:
    real(kind=8) :: xcent,ycent,xoff,yoff
    integer :: iindex,jindex,mq !ioff, joff
    real(kind=8) :: h(4),drytol2,topo

    iindex =  int((xc-xlower)/dx) + 1
    jindex =  int((yc-ylower)/dy) + 1

    xcent  = xlower + (iindex-.5d0)*dx
    ycent  = ylower + (jindex-.5d0)*dy
    xoff   = (xc-xcent)/dx
    yoff   = (yc-ycent)/dy

    drytol2 = 0.1d0 * dry_tolerance

    h(1) = q(1,iindex,jindex)
    h(2) = q(1,iindex+1,jindex)
    h(3) = q(1,iindex,jindex+1)
    h(4) = q(1,iindex+1,jindex+1)

    if ((h(1) < drytol2) .or.  &
        (h(2) < drytol2) .or.  &
        (h(3) < drytol2) .or.  &
        (h(4) < drytol2)) then
        !! One of the cells is dry, so just use value from grid cell
        !! that contains gauge rather than interpolating

        !! icell = int(1.d0 + (xc - xlower) / hx)
        !! jcell = int(1.d0 + (yc - ylower) / hy)
        do mq=1,3
            qvar(mq) = q(mq,iindex,jindex)
        enddo
        !! This is the bottom layer and we should figure out the
        !! topography
        topo = aux(1,iindex,jindex)
    else
        !! Linear interpolation between four cells
        do mq=1,3
            qvar(mq) = (1.d0 - xoff) * (1.d0 - yoff) &
                * q(mq,iindex,jindex)  &
                + xoff*(1.d0 - yoff) * q(mq,iindex+1,jindex)  &
                + (1.d0 - xoff) * yoff * q(mq,iindex,jindex+1)  &
                + xoff * yoff * q(mq,iindex+1,jindex+1)
        enddo
        topo = (1.d0 - xoff) * (1.d0 - yoff)  &
                  * aux(1,iindex,jindex)  &
                  + xoff * (1.d0 - yoff) * aux(1,iindex+1,jindex)  &
                  + (1.d0 - xoff) * yoff * aux(1,iindex,jindex+1)  &
                  + xoff * yoff * aux(1,iindex+1,jindex+1)
    endif

    ! Extract surfaces
    !!eta = qvar(1) + topo  !! Compute this when printing buffers
    avar(1) = topo

    ! Zero out tiny values to prevent later problems reading data,
    ! as done in valout.f
    do mq = 1,3
        if (abs(qvar(mq)) < 1d-99) qvar(mq) = 0.d0
    end do
    !! if (abs(eta) < 1d-99) eta = 0.d0

END SUBROUTINE fc2d_geoclaw_update_gauge
