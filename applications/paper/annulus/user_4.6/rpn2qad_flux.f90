!! This would generally be supplied by the user

SUBROUTINE  rpn2qad_flux(meqn,maux,idir,iface,q,auxvec, & 
                  auxvec_edge,flux)

    IMPLICIT NONE

    INTEGER meqn,maux,idir, iface
    DOUBLE PRECISION q(meqn), flux(meqn)
    DOUBLE PRECISION auxvec(maux), auxvec_edge(maux)

    DOUBLE PRECISION urot
    INTEGER m

!!  # Get cell-centered velocity projected to face 
!!  # 'iface' (in 0,1,2,3) 
    urot = auxvec(2+iface)

!!  #  f(q) = (n dot u)*q
    DO m = 1,meqn
!!      # Don't multiply by edgelength (scaling is done elsewhere)
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

