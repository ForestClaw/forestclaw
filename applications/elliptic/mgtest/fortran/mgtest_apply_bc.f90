subroutine mgtest_fort_apply_bc(blockno, mx, my,mbc,meqn,xlower,ylower, &
    dx,dy,t,intersects_bc,bctype,rhs,g_bc)

    implicit none

    integer blockno, mx,my,mbc,meqn,intersects_bc(0:3),bctype(0:3)
    double precision xlower,ylower,dx,dy,t,g_bc
    double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc)

    !! Dummy argument needed to apply BC
    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision beta(1-mbc:mx+mbc,1-mbc:my+mbc,3)

    integer i,j, iface, idir, i1, ig, j1, jg
    double precision d, h, x, y, g, dx2, dy2
    double precision a,b
    double precision val_beta, grad_beta(2), div_beta_grad_u, flux(0:3)


    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            x = xlower + (i-0.5)*dx
            y = ylower + (i-0.5)*dy
            call mgtest_beta(x,y,val_beta,grad_beta(2))
            beta(i,j,1) = val_beta

            q(i,j) = 0
        end do
    end do

    do i = 2-mbc,mx+mbc
        do j = 2-mbc,my+mbc
            beta(i,j,2) = (beta(i,j,1) + beta(i-1,j,1))/2.d0
            beta(i,j,3) = (beta(i,j,1) + beta(i,j-1,1))/2.d0
        end do
    end do


    do iface = 0,3
        if (intersects_bc(iface) .ne. 0) then
            idir = iface/2   !! direction : 0 or 1

            !! customize this if a,b depend on x
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

            !! Discretize the boundary conditions as : 
            !!
            !!     a*(qI + qg)/2 + b*(qg - qI)/dx = g
            !! 
            !! where gI is the interior q values and qg is             
            !! the ghost cell value.  Solve for qg : 
            !! 
            !!     qg = (g - d1*qI)/d2
            !!
            !!     d1 = (a/2 - g/h)
            !!     d2 = (a/2 + b/h)  !! == d in code below
            !!
            !! Here, we set qg = g/d2
            d = (a/2.d0 + b/h)  
            if (d .eq. 0) then
                write(6,*) 'mgtest_fort_apply_bc : ill-defined BCs'
                stop
            endif

            if (idir .eq. 0) then
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
            elseif (idir .eq. 1) then
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
        endif
    enddo

    !! Ghost cells now all filled in.  Now apply variable coefficent
    !! Laplace operator

    !! This could be done more efficiently
    dx2 = dx*dx
    dy2 = dy*dy
    do i = 1,mx
        do j = 1,my            
            flux(0) = beta(i,j,2)*(q(i,j) - q(i-1,j))/dx
            flux(1) = beta(i+1,j,2)*(q(i+1,j) - q(i,j))/dx
            flux(2) = beta(i,j,3)*(q(i,j) - q(i,j-1))/dy
            flux(3) = beta(i,j+1,3)*(q(i,j+1) - q(i,j))/dy
            div_beta_grad_u = (flux(1) - flux(0))/dx + (flux(3) - flux(2))/dy
            rhs(i,j) = rhs(i,j) - div_beta_grad_u
        end do
    end do 

end subroutine mgtest_fort_apply_bc

