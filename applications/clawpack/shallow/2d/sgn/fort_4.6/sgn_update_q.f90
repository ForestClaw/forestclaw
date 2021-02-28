subroutine sgn_update_q(mx,my,mbc, meqn, mfields,& 
           xlower,ylower,dx,dy,t, dt, maux,aux,q,D)
    implicit none

    integer meqn, mbc, mx, my, maux, mfields
    double precision xlower, ylower, dx, dy, t, dt
    double precision   q(1-mbc:mx+mbc,1-mbc:my+mbc, meqn)
    double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

    !! SGN solution
    double precision D(1-mbc:mx+mbc,1-mbc:my+mbc, mfields)

    DOUBLE PRECISION  grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level

    DOUBLE PRECISION breaking, alpha
    COMMON /common_sgn/ breaking, alpha

    double precision hc, dh(2), db(2), h, deta, src
    integer i,j,m

    do i = 1,my
        do j = 1,mx
            dh(1) = (q(i+1,j,  1) - q(i-1,j,  1))/(2*dx)
            dh(2) = (q(i,  j+1,1) - q(i,  j-1,1))/(2*dy)
            !! write(6,*) i,j, q(i,j,1), D(i,j,1), D(i,j,2)

            !! Assume flat bottom for now
            db(1) = 0  
            db(2) = 0  

            hc = q(i,j,1)
            do m = 1,2
                deta = dh(m) + db(m)
                if (abs(deta) .gt. breaking) then
                    !! If either component in grad eta is greater than breaking, 
                    !! don't update either component
                    continue
                endif
                q(i,j,1+m) = q(i,j,1+m) + dt*hc*(grav/alpha*deta - D(i,j,m))
            end do
        end do
    end do



    return
    end
