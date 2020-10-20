subroutine multigrid_fort_apply_bc_default(blockno, mx, my,mbc,meqn,xlower,ylower, &
    dx,dy,t,intersects_bc,bctype,rhs,g_bc, cons_check, flux_sum)

    implicit none

    external g_bc
    integer blockno, mx,my,mbc,meqn,intersects_bc(0:3),bctype(0:3), cons_check
    double precision xlower,ylower,dx,dy, t, flux_sum(0:3)
    double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc)

    !! Dummy argument needed to apply BC
    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc)

    integer i,j, iface, idir, i1, ig, j1, jg
    double precision d, h, x, y, g, g_bc, dx2, dy2, qlap
    double precision a,b


    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            q(i,j) = 0
        end do
    end do

    do iface = 0,3
        if (intersects_bc(iface) .ne. 0) then
            idir = iface/2   !! direction : 0 or 1

            if (bctype(iface) .eq. 1) then
                a = 1
                b = 0
            elseif (bctype(iface) .eq. 2) then
                a = 0
                b = 1
            endif

            if (idir == 0) then
                h = dx
            else
                h = dy
            endif

            d = (a/2.d0 + b/h)
            if (d .eq. 0) then
                write(6,*) 'multigrid_fort_apply_bc_default : ill-defined BCs'
                stop
            endif

            if (idir == 0) then
                if (iface .eq. 0) then
                    i1 = 1
                    ig = 0
                elseif (iface .eq. 1) then
                    i1 = mx+1
                    ig = mx+1
                endif
                !! location at interface
                x = xlower + (i1 - 1)*dx    

                do j = 1,my
                    y = ylower + (j-0.5)*dy

                    !! inhomogeneity
                    g = g_bc(iface,t,x,y)

                    !! Assume uI == 0
                    q(ig,j) = g/d
                end do
            elseif (idir == 1) then
                if (iface .eq. 2) then
                    j1 = 1
                    jg = 0
                elseif (iface .eq. 3) then
                    j1 = my+1
                    jg = my+1
                endif
                !! location at interface
                y = ylower + (j1 - 1)*dy

                do i = 1,mx
                    x = xlower + (i-0.5)*dx

                    !! inhomogeneity
                    g = g_bc(iface,t,x,y)

                    !! Assume uI == 0
                    q(i,jg) = g/d
                end do
            endif

            !! Ghost cells now all filled in.  Now apply Laplacian

            !! This could be done more efficiently
            dx2 = dx*dx
            dy2 = dy*dy
            do i = 1,mx
                do j = 1,my
                    qlap = (q(i-1,j) - 2*q(i,j) + q(i+1,j))/dx2 + & 
                          (q(i,j-1) - 2*q(i,j) + q(i,j+1))/dy2
                    rhs(i,j) = rhs(i,j) - qlap
                end do
            end do 
        end if 
    end do

end subroutine multigrid_fort_apply_bc_default

