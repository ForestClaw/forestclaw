subroutine sgn_update_q(mx,my,mbc, meqn, mfields,& 
           xlower,ylower,dx,dy,t, dt, maux,aux,q,D)
    implicit none

    integer meqn, mbc, mx, my, maux, mfields
    double precision xlower, ylower, dx, dy, t, dt
    double precision   q(1-mbc:mx+mbc,1-mbc:my+mbc, meqn)
    double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

    !! SGN solution
    double precision D(1-mbc:mx+mbc,1-mbc:my+mbc, mfields)

!!    DOUBLE PRECISION  grav, dry_tolerance, sea_level
!!    COMMON /common_swe/ grav, dry_tolerance, sea_level

    DOUBLE PRECISION breaking, alpha
    COMMON /common_sgn/ breaking, alpha

    DOUBLE PRECISION  grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level


    double precision hc, dh(2), db(2), h, deta, src
    double precision xc,yc,b, grad(2), d2xzb, d2yzb, d2xyzb, dry_tol
    logical wet(-1:1,2)
    integer i,j,m

    dry_tol = dry_tol
    do i = 1,my
        do j = 1,mx
            xc = xlower + (i-0.5)*dx

            dh(1) = (q(i+1,j,  1) - q(i-1,j,  1))/(2*dx)
            dh(2) = (q(i,  j+1,1) - q(i,  j-1,1))/(2*dy)
            !! write(6,*) i,j, q(i,j,1), D(i,j,1), D(i,j,2)

            call sgn_fort_bathy_complete(xc,yc, b, grad, d2xzb, d2yzb, d2xyzb)

            !! Assume flat bottom for now
            db(1) = grad(1)
            db(2) = grad(2) 

            hc = q(i,j,1)
            wet(-1,1) = q(i-1,j,1)  .gt. dry_tol
            wet(0,1)  = q(i,j,1)    .gt. dry_tol
            wet(1,1)  = q(i+1,j,1)  .gt. dry_tol

            wet(-1,2) = q(i,j-1,1) .gt. dry_tol
            wet(0,2)  = q(i,j,1)   .gt. dry_tol
            wet(1,2)  = q(i,j+1,1) .gt. dry_tol
            do m = 1,2
                deta = dh(m) + db(m)
                if (abs(deta) .gt. breaking) then
                    !! If either component in grad eta is greater than breaking, 
                    !! don't update either component
                    continue
                endif
                if (wet(-1,m) .and. wet(0,m) .and. wet(1,m)) then
                    if (m .eq. 1) then
                        q(i,j,1+m) = q(i,j,1+m) + dt*hc*(grav/alpha*deta - D(i,j,m))
                    endif
                endif
            end do
        end do
    end do



    return
    end
