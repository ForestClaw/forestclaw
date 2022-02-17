!! This applies general non-homogeneous boundary conditions by 
!! modifying the right hand side side with non-zero corrections.
!!
!! If cons_check == 0 :  Apply inhomogeneous boundary conditions
!! If cons_check == 1 :  Sum up fluxes around the boundary
!!
!! The general boundary conditions can be written as 
!! 
!!        a(x,y)*q(x,y) + b(x,y)*q_n(x,y) = g(x,y)
!!
!! on the boundary of a square domain.  Below, we only implement the 
!! case a = 1 or b = 1 (but not both).  This can be easily modified following
!! the comments below. 
!!
!! The approach below requires that the term
!!
!!        (a(x,y)/2.d0 + b(x,y)/h)
!!
!! where h = dx or h = dy, is non-zero.
!!

subroutine poisson_fort_apply_bc(blockno, mx, my,mbc,mfields,xlower,ylower, &
    dx,dy,t,intersects_bc,bctype,rhs,g_bc,cons_check,flux_sum)

    implicit none

    integer blockno, mx,my,mbc,mfields,intersects_bc(0:3),bctype(0:3)
    integer cons_check
    double precision xlower,ylower,dx,dy,t,g_bc, flux_sum(mfields)
    double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)

    !! Dummy arrays needed to apply boundary conditions
    double precision qh(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
    double precision beta(1-mbc:mx+mbc,1-mbc:my+mbc,3)

    integer i,j, m, iface, idir, i1, ig, ic, j1, jg, jc
    double precision d, h, x, y, g
    double precision a,b
    double precision val_beta, grad_beta(2), flux(0:3)
    double precision uI, dI

    logical ccheck

    ccheck = cons_check .ne. 0

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            x = xlower + (i-0.5)*dx
            y = ylower + (j-0.5)*dy
            call poisson_fort_beta(x,y,val_beta,grad_beta)
            beta(i,j,1) = val_beta
        end do
    end do

    do i = 2-mbc,mx+mbc
        do j = 2-mbc,my+mbc
            beta(i,j,2) = (beta(i,j,1) + beta(i-1,j,1))/2.d0
            beta(i,j,3) = (beta(i,j,1) + beta(i,j-1,1))/2.d0
        end do
    end do


    do m = 1,mfields

        !! Homogeneous array;  this is overkill - we don't need the entire array.
        do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc
                qh(i,j,m) = 0
            end do
        end do

        do iface = 0,3
            if (intersects_bc(iface) .ne. 0) then
                idir = iface/2   !! direction : 0 or 1

                !! customize this for more general boundary conditions
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
                !!     a*(qI + qg)/2 + b*(qg - qI)/ds = g
                !! 
                !! where gI is the interior q values and qg is             
                !! the ghost cell value.  Solve for qg : 
                !! 
                !!     qg = (g - d1*qI)/d2
                !!
                !!     d1 = (a/2 - b/h)
                !!     d2 = (a/2 + b/h)  !! == d in code below
                !!
                !! Here, we set qg = g/d2
                dI = (a/2.d0 - b/h)  !! Needed for cons check
                d = (a/2.d0 + b/h)  
                if (d .eq. 0) then
                    write(6,*) 'poisson_fort_apply_bc : ill-defined BCs'
                    stop
                endif

                if (idir .eq. 0) then
                    if (iface .eq. 0) then
                        ic = 1  !! cell center
                        i1 = 1
                        ig = 0
                    elseif (iface .eq. 1) then
                        ic = mx  !! cell center
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
                        uI = 0
                        if (ccheck) then
                            uI = rhs(ic,j,m)
                        endif
                        qh(ig,j,m) = (g - dI*uI)/d
                    end do
                elseif (idir .eq. 1) then
                    if (iface .eq. 2) then
                        jc = 1
                        j1 = 1
                        jg = 0
                    elseif (iface .eq. 3) then
                        jc = my
                        j1 = my+1
                        jg = my+1
                    endif
                    !! location at interface
                    y = ylower + (j1 - 1)*dy

                    do i = 1,mx
                        x = xlower + (i-0.5)*dx

                        !! inhomogeneity
                        g = g_bc(iface,t,x,y)

                        uI = 0
                        if (ccheck) then
                            uI = rhs(i,jc,m)
                        endif

                        qh(i,jg,m) = (g - dI*uI)/d
                    end do
                endif
            endif
        end do

        !! Ghost cells now all filled in.  Now apply variable coefficent
        !! Laplace operator

        if (intersects_bc(0) .ne. 0) then
            do j = 1,my
                if (ccheck) then
                    flux(0) = beta(1,j,2)*(rhs(1,j,m) - qh(0,j,m))/dx
                    flux_sum(m) = flux_sum(m) - flux(0)*dy    
                else
                    flux(0) = beta(1,j,2)*(qh(1,j,m) - qh(0,j,m))/dx
                    rhs(1,j,m) = rhs(1,j,m) - (-flux(0)/dx)
                endif
            end do
        endif            
        
        if (intersects_bc(1) .ne. 0) then
            do j = 1,my
                if (ccheck) then
                    flux(1) = beta(mx+1,j,2)*(qh(mx+1,j,m) - rhs(mx,j,m))/dx
                    flux_sum(m) = flux_sum(m) + flux(1)*dy
                else
                    flux(1) = beta(mx+1,j,2)*(qh(mx+1,j,m) - qh(mx,j,m))/dx
                    rhs(mx,j,m) = rhs(mx,j,m) - (flux(1)/dx)
                endif
            end do
        endif

        if (intersects_bc(2) .ne. 0) then
            do i = 1,mx
                if (ccheck) then
                    flux(2) = beta(i,1,3)*(rhs(i,1,m) - qh(i,0,m))/dy
                    flux_sum(m) = flux_sum(m) - flux(2)*dx
                else
                    flux(2) = beta(i,1,3)*(qh(i,1,m) - qh(i,0,m))/dy
                    rhs(i,1,m) = rhs(i,1,m) - (-flux(2)/dy)
                endif
            end do
        endif

        if (intersects_bc(3) .ne. 0) then
            do i = 1,mx
                if (ccheck) then
                    flux(3) = beta(i,my+1,3)*(qh(i,my+1,m) - rhs(i,my,m))/dy           
                    flux_sum(m) = flux_sum(m) + flux(3)*dx
                else
                    flux(3) = beta(i,my+1,3)*(qh(i,my+1,m) - qh(i,my,m))/dy           
                    rhs(i,my,m) = rhs(i,my,m) - (flux(3)/dy)
                endif
            end do
        endif

    end do

end subroutine poisson_fort_apply_bc

