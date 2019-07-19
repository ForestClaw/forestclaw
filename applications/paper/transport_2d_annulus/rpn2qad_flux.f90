!! This would generally be supplied by the user

SUBROUTINE  rpn2qad_flux(meqn,maux,idir,iface,q,auxvec,flux)

    IMPLICIT NONE

    INTEGER meqn,maux,idir, iface
    DOUBLE PRECISION q(meqn), flux(meqn)
    DOUBLE PRECISION auxvec(maux)

    DOUBLE PRECISION urot, g, nv(2), u(2)
    DOUBLE PRECISION m

    g = auxvec(12+iface)  !! Edge length for iface = 0,1,2,3

    u(1) = auxvec(2)
    u(2) = auxvec(3)

    nv(1) = auxvec(4+2*iface)
    nv(2) = auxvec(5+2*iface)

    urot = g*(nv(1)*u(1) + nv(2)*u(2))

    DO m = 1,meqn
        !! # Scaling done here (unlike in ForestClaw)    
        flux(m) = urot*q(m)
    ENDDO

END SUBROUTINE rpn2qad_flux
