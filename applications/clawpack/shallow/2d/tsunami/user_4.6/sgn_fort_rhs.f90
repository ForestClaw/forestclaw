subroutine sgn_fort_rhs(blockno, mbc,mx,my,meqn,mfields, & 
                         xlower,ylower,dx,dy,q,rhs)
    IMPLICIT NONE

    INTEGER mbc,mx,my, mfields, meqn
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

    INTEGER i,j, m, blockno
    DOUBLE PRECISION h, dudx,dudy,dvdx,dvdy, divu
    double precision dc(2),dh(2),db(2),deta, R1, R2

    do j = 1-mbc,mx+mbc
        do i = 1-mbc,my+mbc
            h = q(i,j,1)
            if (h .le. 0) then
                write(6,*) 'fort_rhs : h .le. 0'
                stop
            endif
            u(i,j) = q(i,j,2)/h
            v(i,j) = q(i,j,3)/h
        end do
    end do



    do j = 1,my
        do i = 1,mx
            !! Compute dudx, dvdx,dudy, dvdy
            dudx = (u(i+1,j) - u(i-1,j))/(2*dx)
            dudy = (u(i,j+1) - u(i,j-1))/(2*dy)

            dvdx = (v(i+1,j) - v(i-1,j))/(2*dx)
            dvdy = (v(i,j+1) - v(i,j-1))/(2*dy)

            divu = (dudx + dvdy)
            c(i,j) = -dudx*dvdy + dvdx*dudy + divu**2
        end do
    end do

    do j = 1,my
        do i = 1,mx

            dh(1) = (q(i+1,j,  1) - q(i-1,j,  1))/(2*dx)
            dh(2) = (q(i,  j+1,1) - q(i,  j-1,1))/(2*dy)

            dc(1) = (c(i+1,j) - c(i-1,j))/(2*dx)
            dc(2) = (c(i,j+1) - c(i,j-1))/(2*dy)

            !! Assume flat bottom for now
            db(1) = 0  
            db(2) = 0  

            h = q(i,j,1)
            do m = 1,2
                deta = dh(m) + db(m)

                R1 = -h*(h/3.0*dc(m) + c(i,j)*(dh(m) + db(m)/2.0))                
                R2 = 0    !! Flat bottom assumption

                rhs(i,j,m) = h*(grav/alpha*deta + (-2*R1 + R2))
            end do
            rhs(i,j,2) = 0
        end do
    end do

end subroutine sgn_fort_rhs


