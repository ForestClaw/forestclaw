!! This would generally be supplied by the user

SUBROUTINE  rpn2qad_flux(meqn,maux,idir,iface,q, auxvec,flux)

    IMPLICIT NONE

    INTEGER meqn,maux,idir, m, iface
    DOUBLE PRECISION q(meqn), flux(meqn)
    DOUBLE PRECISION auxvec(maux)

    INTEGER qad_mode
    COMMON /qad_comm/ qad_mode

    DOUBLE PRECISION urot, g, nv(2), u(2)
    INTEGER k

    write(6,*) 'qad_flux is not yet working'
    stop

    !!    Use this version if only 2 lengths (left, bottom) are stored in 
    !!    each cell and aux variables have been identified as "xleft", "yleft", etc.      
    IF (qad_mode .le. 1) THEN
        !! This won't work, but illustrates behavior for original qad
        g = auxvec(6+2*idir) !! Edge length for iface = 0,1,2,3
    ELSE
        !! Use this version if all four edge lengths are available in each cell. 
        g = auxvec(6+iface)  !! Edge length for iface = 0,1,2,3
    ENDIF

    DO m = 1,meqn
        !! # Scaling done here (unlike in ForestClaw)    
        flux(m) = g*urot*q(m)
    ENDDO

END SUBROUTINE rpn2qad_flux
