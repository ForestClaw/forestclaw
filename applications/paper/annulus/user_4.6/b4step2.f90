!> Called before each call to step2.
!! Use to set time-dependent aux arrays or perform other tasks.
!!
!! This default version does nothing. 

subroutine clawpack46_b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,xlower,ylower, &
    dx,dy,t,dt,maux,aux)

    implicit none
    integer maxmx, maxmy, mbc,mx,my,meqn,maux
    double precision  xlower,ylower,dx,dy,t,dt
    double precision  q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    integer i,j
    double precision xc,yc, vel(3)

    !! # Cell-centered velocities : entries (4,5,6)
    DO i = 1-mbc,mx+mbc
        DO j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            CALL annulus_center_velocity(xc,yc,t, vel)

            aux(i,j,2) = vel(1)
            aux(i,j,3) = vel(2)
        end do
    end do


end subroutine clawpack46_b4step2
