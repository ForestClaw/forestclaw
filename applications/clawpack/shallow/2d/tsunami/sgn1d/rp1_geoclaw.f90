!!======================================================================

SUBROUTINE rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,&
               qr,auxl,auxr,fwave,s,amdq,apdq)

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
  INTEGER maxm,meqn,maux,mwaves,mbc,mx

  DOUBLE PRECISION  fwave(meqn, mwaves, 1-mbc:maxm+mbc)
  DOUBLE PRECISION  s(mwaves, 1-mbc:maxm+mbc)
  DOUBLE PRECISION  ql(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION  qr(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION  apdq(meqn,1-mbc:maxm+mbc)
  DOUBLE PRECISION  amdq(meqn,1-mbc:maxm+mbc)
  DOUBLE PRECISION  auxl(maux,1-mbc:maxm+mbc)
  DOUBLE PRECISION  auxr(maux,1-mbc:maxm+mbc)

  DOUBLE PRECISION :: grav, dry_tolerance, sea_level
  COMMON /common_swe/ grav, dry_tolerance, sea_level

  !local only
  INTEGER m,i,mw,maxiter,mu,nv
  DOUBLE PRECISION wall(3)
  DOUBLE PRECISION fw(3,3)
  DOUBLE PRECISION sw(3)

  DOUBLE PRECISION hR,hL,huR,huL,uR,uL,phiR,phiL,pL,pR
  DOUBLE PRECISION bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
  DOUBLE PRECISION s1m,s2m
  DOUBLE PRECISION hstar,hstartest,hstarHLL,sLtest,sRtest
  DOUBLE PRECISION tw,dxdc

  DOUBLE PRECISION  hvR, hvL, vR, vL, rho
  LOGICAL use_simple    

  LOGICAL rare1,rare2

  integer ii_com, jj_com
  common /common_ii/ ii_com, jj_com

  use_simple = .false.

  pL = 0
  pR = 0
  hvR = 0
  hvL = 0
  vR = 0
  vL = 0

  !!g = grav

  rho = -9999999.d0 !! To make sure this doesn't get used. 

  !loop through Riemann problems at each grid cell
  DO i = 2-mbc,mx+mbc
      ii_com = i

      !! -----------------------Initializing-----------------------------------
      !! inform of a bad riemann problem from the start
     IF((qr(1,i-1) .LT. 0.d0) .OR. (ql(1,i) .LT. 0.d0)) THEN
        WRITE(6,201) 'Negative input: hl,hr,i=',i, qr(1,i-1),ql(1,i)
        stop
     ENDIF
201  format(A,I5,2E12.4)          

     !!Initialize Riemann problem for grid interface
     DO mw=1,mwaves
        s(mw,i)=0.d0
        fwave(1,mw,i)=0.d0
        fwave(2,mw,i)=0.d0
        fwave(3,mw,i) = 0
     ENDDO

     !!zero (small) negative values if they exist
     IF (qr(1,i-1) .LT. 0.d0) THEN
        qr(1,i-1)=0.d0
        qr(2,i-1)=0.d0
     ENDIF

     IF (ql(1,i).LT.0.d0) THEN
        ql(1,i)=0.d0
        ql(2,i)=0.d0
     ENDIF

     !!skip problem if in a completely dry area
     IF (qr(1,i-1) <= dry_tolerance .AND. ql(1,i) <= dry_tolerance) THEN
        go to 30
     ENDIF

     !! Riemann problem variables
     hL = qr(1,i-1)
     hR = ql(1,i)
     huL = qr(2,i-1)
     huR = ql(2,i)
     bL = auxr(1,i-1)
     bR = auxl(1,i)

     hvR = 0
     hvL = 0
     vR = 0
     vL = 0

     !!check for wet/dry boundary
     IF (hR .GT. dry_tolerance) THEN
        uR = huR/hR
        phiR = 0.5d0*grav*hR**2 + huR**2/hR
     ELSE
        hR = 0.d0
        huR = 0.d0
        uR = 0.d0
        phiR = 0.d0
     ENDIF

     IF (hL .gt. dry_tolerance) THEN
        uL = huL/hL
        phiL = 0.5d0*grav*hL**2 + huL**2/hL
     ELSE
        hL=0.d0
        huL=0.d0
        uL=0.d0
        phiL = 0.d0
     ENDIF

     wall(1) = 1.d0
     wall(2) = 1.d0
     wall(3) = 1.d0
     IF (hR .LE. dry_tolerance) THEN
        CALL riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m, &
             rare1,rare2,1,dry_tolerance,grav)

        hstartest = MAX(hL,hstar)
        IF (hstartest + bL .LT. bR) THEN
           !!right state should become ghost values that mirror left for wall problem
           !! bR=hstartest+bL
           wall(1)=0.d0
           wall(2)=0.d0
           hR = hL
           huR = -huL
           bR = bL
           phiR = phiL
           uR = -uL
           vL = vR
        ELSEIF (hL+bL.LT.bR) THEN
           bR = hL + bL
        ENDIF
     ELSEIF (hL .LE. dry_tolerance) THEN ! right surface is lower than left topo
        CALL riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m, &
             rare1,rare2,1,dry_tolerance,grav)
        hstartest = MAX(hR,hstar)
        IF (hstartest + bR .LT. bL) THEN
           !!left state should become ghost values that mirror right
           !! bL=hstartest+bR
           wall(1) = 0.d0
           wall(2) = 0.d0
           hL = hR
           huL = -huR
           bL = bR
           phiL = phiR
           uL = -uR
           vR = vL
        ELSEIF (hR + bR .LT. bL) THEN
           bL = hR + bR
        ENDIF
     ENDIF

     !!determine wave speeds
     sL = uL - SQRT(grav*hL) ! 1 wave speed of left state
     sR = uR + SQRT(grav*hR) ! 2 wave speed of right state

     uhat = (SQRT(grav*hL)*uL + SQRT(grav*hR)*uR)/(SQRT(grav*hR) + SQRT(grav*hL)) ! Roe average
     chat = SQRT(grav*0.5d0*(hR + hL)) ! Roe average
     sRoe1 = uhat - chat ! Roe wave speed 1 wave
     sRoe2 = uhat + chat ! Roe wave speed 2 wave

     sE1 = MIN(sL,sRoe1) ! Eindfeldt speed 1 wave
     sE2 = MAX(sR,sRoe2) ! Eindfeldt speed 2 wave

     !!--------------------end initializing...finally----------

     !!solve Riemann problem.

     maxiter = 1

     if (use_simple) then                 
        CALL  simple_riemann(hR,uR,vr, hL,uL,vl, uhat,chat,bL, bR, &
                             phiR,phiL,sw,fw)
    else
        jj_com = 1
        CALL riemann_aug_JCP(maxiter,3,3,hL,hR,huL, huR, & 
                             hvL,hvR,bL,bR,uL,uR,vL,vR, phiL, phiR, &
                             pL,pR,sE1,sE2,dry_tolerance,grav,rho,sw,fw)
    endif

     !! eliminate ghost fluxes for wall
     DO mw = 1,mwaves
        sw(mw)   = sw(mw)*wall(mw)
        fw(1,mw) = fw(1,mw)*wall(mw)
        fw(2,mw) = fw(2,mw)*wall(mw)
        fw(3,mw) = fw(3,mw)*wall(mw)
     ENDDO

     DO mw = 1,mwaves
        s(mw,i) = sw(mw)
        fwave(1,mw,i) = fw(1,mw)
        fwave(2,mw,i) = fw(2,mw)
        fwave(3,mw,i) = fw(3,mw)
     ENDDO

30   CONTINUE
  ENDDO

  !!===============================================================================


  !!============= compute fluctuations=============================================
  amdq(1:3,:) = 0.d0
  apdq(1:3,:) = 0.d0
  DO i = 2-mbc,mx+mbc
     DO  mw = 1,3
        IF (s(mw,i) < 0.d0) THEN
           amdq(1:3,i) = amdq(1:3,i) + fwave(1:3,mw,i)
        ELSE IF (s(mw,i) > 0.d0) THEN
           apdq(1:3,i)  = apdq(1:3,i) + fwave(1:3,mw,i)
        ELSE
           amdq(1:3,i) = amdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
           apdq(1:3,i) = apdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
        ENDIF
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE rp1


subroutine simple_riemann(hr,ur,vr, hl,ul,vl, uhat,chat,bl, br, &
                 phir,phil,s,fwave)
    implicit none

    double precision hr,ur,vr, hl,ul,vl, uhat, chat, phir, &
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
        
    !! # Right eigenvectors (column)
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

end


