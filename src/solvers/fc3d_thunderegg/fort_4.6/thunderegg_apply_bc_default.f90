subroutine fc3d_thunderegg_fort_apply_bc_default(blockno, mx, my, mz, mbc,meqn,xlower,ylower,zlower, &
    dx,dy,dz,t,intersects_bc,bctype,rhs,g_bc, cons_check, flux_sum)

    implicit none

    external g_bc
    integer blockno, mx,my,mz,mbc,meqn,intersects_bc(0:5),bctype(0:5), cons_check
    double precision xlower,ylower,zlower,dx,dy,dz, t, flux_sum(0:5)
    double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

    !! Dummy argument needed to apply BC
    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

    integer i,j,k, iface, idir, i1, ig, j1, jg, k1, kg
    double precision d, h, x, y, z, g, g_bc, dx2, dy2, dz2, qlap
    double precision a,b


    do k = 1-mbc,mz+mbc
        do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
                q(i,j,k) = 0
            end do
        end do
    end do

    do iface = 0,5
        if (intersects_bc(iface) .ne. 0) then
            idir = iface/2   !! direction : 0, 1, or 2

            if (bctype(iface) .eq. 1) then
                a = 1
                b = 0
            elseif (bctype(iface) .eq. 2) then
                a = 0
                b = 1
            else
                a = 0
                b = 0
            endif

            if (idir == 0) then
                h = dx
            else if (idir == 1) then
                h = dy
            else
                h = dz
            endif

            d = (a/2.d0 + b/h)
            if (d .eq. 0) then
                write(6,*) 'fc3d_thunderegg_fort_apply_bc_default : ill-defined BCs'
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

                do k = 1,mz
                    z = zlower + (k-0.5)*dz
                    do j = 1,my
                        y = ylower + (j-0.5)*dy

                        !! inhomogeneity
                        g = g_bc(iface,t,x,y,z)

                        !! Assume uI == 0
                        q(ig,j,k) = g/d
                    end do
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

                do k = 1,mz
                    z = zlower + (k-0.5)*dz
                    do i = 1,mx
                        x = xlower + (i-0.5)*dx

                        !! inhomogeneity
                        g = g_bc(iface,t,x,y,z)

                        !! Assume uI == 0
                        q(i,jg,k) = g/d
                    end do
                end do
            elseif (idir == 2) then
                if (iface .eq. 4) then
                    k1 = 1
                    kg = 0
                elseif (iface .eq. 5) then
                    k1 = mz+1
                    kg = mz+1
                endif
                !! location at interface
                z = zlower + (k1 - 1)*dz

                do j = 1,my
                    y = ylower + (j-0.5)*dy
                    do i = 1,mx
                        x = xlower + (i-0.5)*dx

                        !! inhomogeneity
                        g = g_bc(iface,t,x,y,z)

                        !! Assume uI == 0
                        q(i,j,kg) = g/d
                    end do
                end do
            endif

            !! Ghost cells now all filled in.  Now apply Laplacian

            !! This could be done more efficiently
            dx2 = dx*dx
            dy2 = dy*dy
            dz2 = dz*dz
            do k = 1,mz
                do j = 1,my
                    do i = 1,mx
                        qlap = (q(i-1,j,k) - 2*q(i,j,k) + q(i+1,j,k))/dx2 + & 
                               (q(i,j-1,k) - 2*q(i,j,k) + q(i,j+1,k))/dy2 + &
                               (q(i,j,k-1) - 2*q(i,j,k) + q(i,j,k+1))/dz2
                        rhs(i,j,k) = rhs(i,j,k) - qlap
                    end do
                end do
            end do 
        end if 
    end do

end subroutine fc3d_thunderegg_fort_apply_bc_default

