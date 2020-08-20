!! ============================================
SUBROUTINE fc2d_geoclaw_b4step2(mbc,mx,my,meqn,q,xlower,ylower, &
    dx,dy,t,dt,maux,aux)
    !! ============================================
    !!
    !! # called before each call to step
    !! # use to set time-dependent aux arrays or perform other tasks.
    !!
    !! This particular routine sets negative values of q(1,i,j) to zero,
    !! as well as the corresponding q(m,i,j) for m=1,meqn.
    !! This is for problems where q(1,i,j) is a depth.
    !! This should occur only because of rounding error.
    !!
    !! Also calls movetopo if topography might be moving.

    USE geoclaw_module, ONLY: dry_tolerance
    !!USE geoclaw_module, ONLY: g => grav
    USE topo_module, ONLY: num_dtopo  !!,topotime
    USE topo_module, ONLY: aux_finalized
    !!USE topo_module, ONLY: xlowdtopo,xhidtopo,ylowdtopo,yhidtopo
    USE topo_module, ONLY: t0dtopo, tfdtopo

    !!USE amr_module, ONLY: xlowdomain => xlower
    !!USE amr_module, ONLY: ylowdomain => ylower
    !!USE amr_module, ONLY: xhidomain => xupper
    !!USE amr_module, ONLY: yhidomain => yupper
    !!USE amr_module, ONLY: xperdom,yperdom,spheredom,NEEDS_TO_BE_SET
    USE amr_module, ONLY: NEEDS_TO_BE_SET

    USE storm_module, ONLY: set_storm_fields

    IMPLICIT NONE

    !! Subroutine arguments
    INTEGER, INTENT(in) :: meqn
    INTEGER, INTENT(inout) :: mbc,mx,my,maux
    REAL(kind=8), INTENT(inout) :: xlower, ylower, dx, dy, t, dt
    REAL(kind=8), INTENT(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    REAL(kind=8), INTENT(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    !! Local storage
    INTEGER :: index,i,j,k,dummy
    REAL(kind=8) :: h,u,v

    INTEGER :: is_ghost, mint, nghost
    REAL(KIND=8) :: tmin, tmax
    LOGICAL :: dtint1, dtint2

    !! Check for NaNs in the solution
    CALL check4nans(meqn,mbc,mx,my,q,t,1)

    !! check for h < 0 and reset to zero
    !! check for h < dry tolerance
    !! set hu = hv = 0 in all these cells
    FORALL(i=1-mbc:mx+mbc, j=1-mbc:my+mbc,q(1,i,j) < dry_tolerance)
        q(1,i,j) = MAX(q(1,i,j),0.d0)
        q(2:3,i,j) = 0.d0
    END FORALL

    !! Get interval bracketing largest dtopo interval.  
    tmin = 1.d99
    tmax = 0
    do i = 1,num_dtopo
        if (t0dtopo(i) < tmin) then
            tmin = t0dtopo(i)
        endif
        if (tfdtopo(i) > tmax) then
            tmax = tfdtopo(i)
        endif
    enddo

    !! dtopo time is an interval [t0,tf]
    dtint1 = (tmin .le. t .and. t .le. tmax) 

    !! dtopo time is instanteous;; [t0, \infty]
    dtint2 = (tmin .eq. tmax) .and. t .ge. tmax

    if (dtint1 .or. dtint2) then
        !! topo arrays might have been updated by dtopo more recently than
        !! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting

        is_ghost = 0
        nghost = mbc    !! won't be used, if is_ghost = 0
        mint = 2*mbc    !! not used
        CALL fc2d_geoclaw_setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux,is_ghost,nghost,mint)
    ENDIF

    !! Set wind and pressure aux variables for this grid
    CALL set_storm_fields(maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)

END SUBROUTINE fc2d_geoclaw_b4step2
