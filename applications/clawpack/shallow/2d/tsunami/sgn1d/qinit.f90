SUBROUTINE qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
    IMPLICIT NONE

    INTEGER meqn,mbc,mx,maux
    DOUBLE PRECISION xlower,dx

    DOUBLE PRECISION q(meqn,1-mbc:mx+mbc)
    DOUBLE PRECISION aux(maux,1-mbc:mx+mbc)

    DOUBLE PRECISION  grav
    COMMON /common_swe/ grav

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    double precision :: a,b,h0
    common /common_initheight/ a,b,h0

    DOUBLE PRECISION xc,hr, hl, ur
    INTEGER i


!!      
    hl = h0 + a
    hr = h0
    ur = 0.d0

    do i = 1,mx
        xc = xlower + (i-0.5d0)*dx
        if (abs(xc) .lt. b) then
            q(1,i) = hl
        else
            q(1,i) = hr
        endif
        q(2,i) = 0
    end do

    return
end subroutine qinit

