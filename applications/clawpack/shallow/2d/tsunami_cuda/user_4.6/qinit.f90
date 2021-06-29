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

    double precision :: a,b,h0
    common /common_initheight/ a,b,h0

    DOUBLE PRECISION xc,yc,hl, hr, ur, ybar
    INTEGER i,j

!!      
    hl = h0 + a
    hr = h0
    ur = 0.d0

    ybar = ylower + my*dy/2.d0


    do i = 1-mbc,mx+mbc
        xc = xlower + (i-0.5)*dx  
        do j = 1-mbc,my+mbc
            yc = ylower + (j-0.5)*dy
            if (abs(xc) .lt. b) then
                q(i,j,1) = h0 + a*exp(-20*(yc-ybar)**2)
            else
                q(i,j,1) = hr
            endif
            q(i,j,2) = 0
            q(i,j,3) = 0
        end do
    end do

    return
end subroutine clawpack46_qinit

