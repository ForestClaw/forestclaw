subroutine sgn_fort_rhs(blockno, mbc,mx,my,meqn,mfields, & 
                         xlower,ylower,dx,dy,q,rhs)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn, blockno
    DOUBLE PRECISION xlower,ylower,dx,dy
    DOUBLE PRECISION rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)    
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    DOUBLE PRECISION breaking, alpha
    COMMON /common_sgn/ breaking, alpha

    DOUBLE PRECISION  grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level

    !! Dummy variables
    DOUBLE PRECISION u(1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION v(1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION c(1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION d(1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER i,j, m
    DOUBLE PRECISION h, dudx,dudy,dvdx,dvdy, divu
    double precision dc(2),dh(2),db(2),deta, R1, R2, r0_norm

    double precision xc,yc, grad(2), d2xzb, d2yzb, d2xyzb
    double precision dxzb, dyzb, dd(2), b

    rhs = 0
    do j = 1-mbc,mx+mbc
        do i = 1-mbc,my+mbc
            h = q(i,j,1)
            if (h .le. 0) then
                !! write(6,*) 'fort_rhs : h .le. 0'
                !! stop
                u(i,j) = 0
                v(i,j) = 0
            else
                u(i,j) = q(i,j,2)/h
                v(i,j) = q(i,j,3)/h
            endif
        end do
    end do

    do j = 0,my+1
        do i = 0,mx+1
            !! Compute dudx, dvdx,dudy, dvdy
            dudx = (u(i+1,j) - u(i-1,j))/(2*dx)
            dudy = (u(i,j+1) - u(i,j-1))/(2*dy)

            dvdx = (v(i+1,j) - v(i-1,j))/(2*dx)
            dvdy = (v(i,j+1) - v(i,j-1))/(2*dy)

            call sgn_fort_bathy_complete(xc,yc, b,grad, d2xzb, d2yzb, d2xyzb)

            dxzb = grad(1)
            dyzb = grad(2)

            divu = (dudx + dvdy)
            c(i,j) = -dudx*dvdy + dvdx*dudy + divu**2
            d(i,j) = u(i,j)**2*d2xzb + v(i,j)**2*d2yzb + & 
                     2*u(i,j)*v(i,j)*d2xyzb
        end do
    end do

    r0_norm = 0
    do j = 1,my
        do i = 1,mx

            dh(1) = (q(i+1,j,  1) - q(i-1,j,  1))/(2*dx)
            dh(2) = (q(i,  j+1,1) - q(i,  j-1,1))/(2*dy)

            dc(1) = (c(i+1,j) - c(i-1,j))/(2*dx)
            dc(2) = (c(i,j+1) - c(i,j-1))/(2*dy)

            dd(1) = (d(i+1,j) - d(i-1,j))/(2*dx)
            dd(2) = (d(i,j+1) - d(i,j-1))/(2*dy)

            call sgn_fort_bathy_complete(xc,yc, b, grad, d2xzb, d2yzb, d2xyzb)

            !! Assume flat bottom for now
            db(1) = grad(1)
            db(2) = grad(2) 

            h = q(i,j,1)
            do m = 1,2
                deta = dh(m) + db(m)

                R1 = -h*(h/3.*dc(m) + c(i,j)*(dh(m) + db(m)/2.0))                
                R2 =    (h/2.*dd(m) + d(i,j)*(dh(m) + db(m)))    !! Flat bottom assumption

                if (h .gt. 0) then
                    rhs(i,j,m) = h*(grav/alpha*deta + (-2*R1 + R2))
                else
                    rhs(i,j,m) = 0;
                endif
            end do
            if (.true.) then
                !! For pseudo-1d
                rhs(i,j,2) = 0
            endif
        end do
    end do


end subroutine sgn_fort_rhs


