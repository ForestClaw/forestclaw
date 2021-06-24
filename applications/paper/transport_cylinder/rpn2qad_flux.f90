!! This would generally be supplied by the user

SUBROUTINE  rpn2qad_flux(meqn,maux,idir,iface,q, auxvec,flux)

    IMPLICIT NONE

    INTEGER meqn,maux,idir, m, iface
    DOUBLE PRECISION q(meqn), flux(meqn)
    DOUBLE PRECISION auxvec(maux)

    DOUBLE PRECISION urot, g
    INTEGER k

!!  # Cell-centered velocity has been projected to normals on each face.
    urot = auxvec(2+iface)  

    g = auxvec(6+iface)  !! Edge length for iface = 0,1,2,3

    DO m = 1,meqn
        !! # Scaling done here (unlike in ForestClaw)    
        flux(m) = g*urot*q(m)
    ENDDO

END SUBROUTINE rpn2qad_flux
