    SUBROUTINE clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my, & 
    xlower,ylower,dx,dy,q,maux,aux)

    IMPLICIT NONE

    INTEGER maxmx, maxmy, meqn,mbc,mx,my, maux
    DOUBLE PRECISION xlower,ylower, dx,dy

    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    DOUBLE PRECISION aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    DOUBLE PRECISION  grav
    COMMON /common_swe/ grav

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION  grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level


    DOUBLE PRECISION xc,yc,hl, hr, ur, ybar
    INTEGER i,j

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            eta = sea_level
            q(i,j,1) = max(sea_level, eta - aux(i,j,1))
            q(i,j,2) = 0
            q(i,j,3) = 0
        end do
    end do

    return
end subroutine clawpack46_qinit

