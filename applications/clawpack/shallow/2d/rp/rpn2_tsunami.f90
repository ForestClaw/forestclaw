!!======================================================================

SUBROUTINE clawpack46_rpn2(ixy,maxm,meqn,mwaves,mbc,mx, & 
           ql, qr,auxl,auxr,fwave,s,amdq,apdq)

!!======================================================================
!!
!! Solves normal Riemann problems for the 2D SHALLOW WATER equations
!!     with topography:
!!     #        h_t + (hu)_x + (hv)_y = 0                           #
!!     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
!!     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #
!!
!! On input, ql contains the state vector at the left edge of each cell
!!     qr contains the state vector at the right edge of each cell
!!
!! This data is along a slice in the x-direction if ixy=1
!!     or the y-direction if ixy=2.
!!
!!  Note that the i'th Riemann problem has left state qr(i-1,:)
!!     and right state ql(i,:)
!!  From the basic clawpack routines, this routine is called with
!!     ql = qr
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water equations.            !
!                                                                           !
!       It allows the user to easily select a Riemann solver in             !
!       riemannsolvers_geo.f. this routine initializes all the variables    !
!       for the shallow water equations, accounting for wet dry boundary    !
!       dry cells, wave speeds etc.                                         !
!                                                                           !
!           David George, Vancouver WA, Feb. 2009                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!  use geoclaw_module, only: g => grav, dry_tolerance => dry_tolerance, rho
!!  use geoclaw_module, only: earth_radius, deg2rad
!!  use amr_module, only: mcapa
!!
!!  use storm_module, only: pressure_forcing, pressure_index


    IMPLICIT NONE

    !input
    INTEGER maxm,meqn,mwaves,mbc,mx, ixy

    DOUBLE PRECISION  fwave(1-mbc:maxm+mbc,meqn, mwaves)
    DOUBLE PRECISION      s(1-mbc:maxm+mbc,mwaves)
    DOUBLE PRECISION     ql(1-mbc:maxm+mbc,meqn)
    DOUBLE PRECISION     qr(1-mbc:maxm+mbc,meqn)
    DOUBLE PRECISION   apdq(1-mbc:maxm+mbc,meqn)
    DOUBLE PRECISION   amdq(1-mbc:maxm+mbc,meqn)
    DOUBLE PRECISION   auxl(1-mbc:maxm+mbc,*)
    DOUBLE PRECISION   auxr(1-mbc:maxm+mbc,*)

    DOUBLE PRECISION :: grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level

    !local only
    INTEGER m,i,mw,mu,mv
    DOUBLE PRECISION wall(3)
    DOUBLE PRECISION fw(3,3)
    DOUBLE PRECISION sw(3)

    DOUBLE PRECISION hR,hL,huR,huL,uR,uL,phiR,phiL,pL,pR
    DOUBLE PRECISION bR,bL,sL,sR,uhat,chat

    DOUBLE PRECISION  hvR, hvL, vR, vL
    double precision z

    integer icom, jcom
    double precision dtcom, dxcom, dycom, tcom
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    double precision szm(3), szp(3)
    integer mq

    if (ixy .eq. 1) then
        mu=2
        mv=3
    else
        mu=3
        mv=2
    endif

    !loop through Riemann problems at each grid cell
    DO i = 2-mbc,mx+mbc

        !! Get Riemann variables Riemann problem variables
        hL  = qr(i-1,1)
        hR  = ql(i,1)

        huL = qr(i-1,mu)
        huR = ql(i,mu)

        bL  = auxr(i-1,1)
        bR  = auxl(i,1)

        hvL = qr(i-1,mv) 
        hvR = ql(i,mv)        

        uR = huR/hR
        vR = hvR/hR
        phiR = 0.5d0*grav*hR**2 + huR**2/hR

        uL = huL/hL
        vL = hvL/hL
        phiL = 0.5d0*grav*hL**2 + huL**2/hL

        !! # Determine wave speeds
        sL = uL - SQRT(grav*hL) !! 1 wave speed of left state
        sR = uR + SQRT(grav*hR) !! 2 wave speed of right state

        uhat = (SQRT(grav*hL)*uL + SQRT(grav*hR)*uR)/(SQRT(grav*hR) + SQRT(grav*hL)) 
        chat = SQRT(grav*0.5d0*(hR + hL)) 

        CALL  simple_riemann(hR,uR,vR, hL,uL,vl, uhat,chat,bL, bR, &
                             phiR,phiL,sw,fw)

        DO mw = 1,mwaves
            s(i,mw) = sw(mw)
            fwave(i,1,mw) = fw(1,mw)
            fwave(i,2,mw) = fw(mu,mw)
            fwave(i,3,mw) = fw(mv,mw)
        ENDDO


30  CONTINUE
    ENDDO

    !!===============================================================================


    !!============= compute fluctuations=============================================
!!    amdq(:,1:3) = 0.d0
!!    apdq(:,1:3) = 0.d0
!!    DO i = 2-mbc,mx+mbc
!!        DO  mw = 1,3
!!!!            IF (s(i,mw) < 0.d0) THEN
!!!!                amdq(i,1:3) = amdq(i,1:3) + fwave(i,1:3,mw)
!!!!            ELSE IF (s(i,mw) > 0.d0) THEN
!!!!                apdq(i,1:3)  = apdq(i,1:3) + fwave(i,1:3,mw)
!!!!            ELSE
!!!!                amdq(i,1:3) = amdq(i,1:3) + 0.5d0 * fwave(i,1:3,mw)
!!!!                apdq(i,1:3) = apdq(i,1:3) + 0.5d0 * fwave(i,1:3,mw)
!!!!            ENDIF
!!            !! if s == 0, returns 1
!!            z = sign(1.d0,s(i,mw))
!!            amdq(i,1:3) = amdq(i,1:3) + (1-z)/2.0*fwave(i,1:3,mw)
!!            apdq(i,1:3) = apdq(i,1:3) + (z+1)/2.0*fwave(i,1:3,mw)
!!        ENDDO
!!    ENDDO

    do i = 2-mbc,mx+mbc
        do mw = 1,mwaves
            z = sign(1.d0,s(i,mw))
            szm(mw) = (1-z)/2
            szp(mw) = (1+z)/2
        end do

        do mq = 1,meqn
            amdq(i,mq) =              szm(1)*fwave(i,mq,1)
            amdq(i,mq) = amdq(i,mq) + szm(2)*fwave(i,mq,2)
            amdq(i,mq) = amdq(i,mq) + szm(3)*fwave(i,mq,3)

            apdq(i,mq) =              szp(1)*fwave(i,mq,1)
            apdq(i,mq) = apdq(i,mq) + szp(2)*fwave(i,mq,2)
            apdq(i,mq) = apdq(i,mq) + szp(3)*fwave(i,mq,3)
        enddo 
    end do


    RETURN
END SUBROUTINE clawpack46_rpn2


SUBROUTINE simple_riemann(hr,ur,vr, hl,ul,vl, uhat,chat,bl, br, &
                 phir,phil,s,fwave)
    IMPLICIT NONE

    DOUBLE PRECISION hr,ur,vr, hl,ul,vl, uhat, chat, phir, &
               phil,s(3), fwave(3,3), bl, bR

    double precision fluxdiff(3),beta(3), hbar

    DOUBLE PRECISION :: grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level

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
        
    !! Wave strengths
    beta(1) = -(fluxdiff(2) - s(3) * fluxdiff(1)) / (s(3) - s(1))
    beta(3) =  (fluxdiff(2) - s(1) * fluxdiff(1)) / (s(3) - s(1))
    beta(2) =   fluxdiff(3) - beta(1)*vl - beta(3)*vr

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


