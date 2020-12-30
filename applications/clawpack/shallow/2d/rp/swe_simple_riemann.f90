SUBROUTINE simple_riemann(hr,ur,vr, hl,ul,vl, uhat,chat,bl, br, &
                 phir,phil,s,fwave)
    IMPLICIT NONE

    DOUBLE PRECISION hr,ur,vr, hl,ul,vl, uhat, chat, phir, &
               phil,s(3), fwave(3,3), bl, bR

    double precision fluxdiff(3),R(3,3), beta(3), hbar

    DOUBLE PRECISION :: grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level

    INTEGER ii_com, jj_com
    COMMON /common_ii/ ii_com, jj_com

    fwave = 0
    s = 0

    hbar = 0.5 * (hr + hl)

    !! # Flux differences
    fluxdiff(1) = (hr * ur) - (hl * ul)
    fluxdiff(2) = phir - phil + grav * hbar * (br - bl)
    fluxdiff(3) = hr * ur * vr - hl * ul * vl

    !! # Wave speeds
    s(1) = MIN(ul - SQRT(grav * hl), uhat - chat)
    s(3) = MAX(ur + SQRT(grav * hr), uhat + chat)
    s(2) = 0.5d0 * (s(1) + s(3))
        
    !! # Right eigenvectors (column) (don't seem to be used)
    R(1,1) = 1.d0
    R(2,1) = s(1)
    R(3,1) = vl
        
    R(1, 2) = 0.d0
    R(2, 2) = 0.0
    R(3, 2) = 1.0

    R(1,3) = 1.d0
    R(2,3) = s(3)
    R(3,3) = vr
    
    !! Wave strengths
    beta(1) = (s(3) * fluxdiff(1) - fluxdiff(2)) / (s(3) - s(1))
    beta(3) = (fluxdiff(2) - s(1) * fluxdiff(1)) / (s(3) - s(1))
    beta(2) = fluxdiff(3) - beta(1)*vl - beta(3)*vr

    !! # Flux waves = beta*R
    fwave(1,1) = beta(1)
    fwave(2,1) = beta(1)*s(1)
    fwave(3,1) = beta(1)*vl

    fwave(1,2) = 0
    fwave(2,2) = 0
    fwave(3,2) = beta(2)


    fwave(1,3) = beta(3)
    fwave(2,3) = beta(3)*s(3)
    fwave(3,3) = beta(3)*vr

END subroutine simple_riemann






