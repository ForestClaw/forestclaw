!! This would generally be supplied by the user

SUBROUTINE  rpn2qad_flux(meqn,maux,idir,iface,q,auxvec, & 
                  auxvec_edge,flux)

    IMPLICIT NONE

    INTEGER meqn,maux,idir, iface
    DOUBLE PRECISION q(meqn), flux(meqn)
    DOUBLE PRECISION auxvec(maux), auxvec_edge(maux)

    DOUBLE PRECISION urot, nv(2), u(2)
    INTEGER m

    u(1) = auxvec(2)
    u(2) = auxvec(3)

    nv(1) = auxvec(4+2*iface)
    nv(2) = auxvec(5+2*iface)

    urot = nv(1)*u(1) + nv(2)*u(2)

    DO m = 1,meqn
        !! Don't scale the flux here 
        flux(m) = urot*q(m)
    ENDDO

END SUBROUTINE rpn2qad_flux

SUBROUTINE  rpn2qad_zero(meqn,maux, idir, iface, q, & 
                         auxvec, auxvec_edge,flux)
    IMPLICIT NONE

    INTEGER meqn,maux,idir, iface
    DOUBLE PRECISION q(meqn), flux(meqn)
    DOUBLE PRECISION auxvec(maux), auxvec_edge(maux)
    INTEGER m

    DO m = 1,meqn
        flux(m) = 0
    END DO

END  SUBROUTINE RPN2QAD_ZERO

