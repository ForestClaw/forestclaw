subroutine clawpack46_qinit(meqn,mbc, mx,my,mz, & 
            xlower,ylower,zlower, dx,dy,dz, q,maux,aux)
    use setprob_mod, only : qin, qout, init_choice, mapping
    implicit none

    integer meqn, mbc, mx, my, mz, maux
    double precision xlower, ylower,zlower,dx,dy,dz

    double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc, meqn)
    double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc, maux)

    integer i,j, k, m, nfix
    double precision xc, yc,zc,xlow,ylow,zlow, win

    integer blockno, fc3d_clawpack46_get_block

    double precision C

    blockno = fc3d_clawpack46_get_block()

    nfix = mx
    if (init_choice .eq. 0) then
        if (mx .ne.  nfix .or. my  .ne. nfix .or. mz .ne. nfix) then
            write(6,*) 'Set mx = my = mz for debugging'
            stop
        endif
        do m = 1,meqn
            q(:,:,:,m) = qout(m)
        enddo
        !! Set one cell to qin values (used in debugging)
        m = nfix/2
        q(m:m+1,m:m+1,m:m+1,1) = qin(1)
        q(m:m+1,m:m+1,m:m+1,5) = qin(5)
    else
        do k = 1-mbc,mz+mbc
            zc = zlower + (k-0.5)*dz
            zlow = zlower + (k-1)*dz
            do  i = 1-mbc,mx+mbc
                xc = xlower + (i-0.5)*dx
                xlow = xlower + (i-1)*dx
                do  j = 1-mbc,my+mbc
                    yc = ylower + (j-0.5)*dy
                    ylow = ylower + (j-1)*dy

                    call cellave3(blockno,xlow,ylow,zlow,dx,dy,dz,win)

                    q(i,j,k,1) = qin(1)*win + (1-win)*qout(1)
                    q(i,j,k,2) = 0.d0
                    q(i,j,k,3) = 0.d0
                    q(i,j,k,4) = 0.d0
                    q(i,j,k,5) = qin(5)*win + (1-win)*qout(5)
                end do
            end do
        end do
    endif

    return

    C = 0
    do i = 1,mx
        do j = 1,my
            do k = 1,mz
                if (mapping .gt. 0) then
                    C = C + q(i,j,k,1)*dx*dy*dz*aux(1,1,1,1)
                else
                    C = C + q(i,j,k,1)*dx*dy*dz
                endif
            end do
        end do
    end do
    write(6,*) 'Total mass = ',C
    write(6,*) qin(1), qout(1)
    if (mapping .gt. 0) then
        write(6,*) aux(1,1,1, 1), dx, dy, dz
    endif
    stop

    return
end subroutine clawpack46_qinit
