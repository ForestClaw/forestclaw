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
    !! If topography has moved, then this also updates the aux arrays with new topography.
    !! 

    USE geoclaw_module, ONLY: dry_tolerance
    USE topo_module, ONLY: num_dtopo  !!,topotime
    USE topo_module, ONLY: aux_finalized
    USE topo_module, ONLY: t0dtopo, tfdtopo

    USE amr_module, ONLY: NEEDS_TO_BE_SET

    !!USE storm_module, ONLY: set_storm_fields

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
    LOGICAL :: fc2d_geoclaw_check_dtopotime, t_in_dtopo_interval

    double precision :: tau


    !! Check for NaNs in the solution
    CALL check4nans(meqn,mbc,mx,my,q,t,1)

    !! check for h < 0 and reset to zero
    !! check for h < dry tolerance
    !! set hu = hv = 0 in all these cells
    FORALL(i=1-mbc:mx+mbc, j=1-mbc:my+mbc,q(1,i,j) < dry_tolerance)
        q(1,i,j) = MAX(q(1,i,j),0.d0)
        q(2:3,i,j) = 0.d0
    END FORALL

    !! Determine if time is in dtopo interval.  If so, we need to update the 
    !! aux array.  
    t_in_dtopo_interval = fc2d_geoclaw_check_dtopotime(t, tau)


    if (t_in_dtopo_interval) then
        ! topo arrays might have been updated by dtopo more recently than
        ! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting
        is_ghost = 0
        nghost = mbc    !! won't be used, if is_ghost = 0
        mint = 2*mbc    !! not used
        CALL fc2d_geoclaw_setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux,is_ghost,nghost,mint)
    endif

    !! Set wind and pressure aux variables for this grid
    !!CALL set_storm_fields(maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)

END SUBROUTINE fc2d_geoclaw_b4step2
    


DOUBLE PRECISION FUNCTION fc2d_geoclaw_get_dt_max_dtopo()
    USE topo_module, ONLY: dt_max_dtopo
    IMPLICIT NONE

    fc2d_geoclaw_get_dt_max_dtopo = dt_max_dtopo

END FUNCTION fc2d_geoclaw_get_dt_max_dtopo



SUBROUTINE fc2d_geoclaw_get_dtopo_interval(tmin, tmax)
    USE topo_module, ONLY: t0dtopo, tfdtopo, num_dtopo

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(OUT) ::  tmin, tmax
    integer :: i

    tmin = 1.d99
    tmax = 0
    do i = 1,num_dtopo
        if (t0dtopo(i) < tmin) then
            tmin = t0dtopo(i)
        endif
        if (tfdtopo(i) >= tmax) then
            tmax = tfdtopo(i)
        endif
    enddo

END SUBROUTINE fc2d_geoclaw_get_dtopo_interval


LOGICAL FUNCTION fc2d_geoclaw_check_dtopotime(t, tau)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) ::  t
    DOUBLE PRECISION, INTENT(out) :: tau

    INTEGER :: i    

    DOUBLE PRECISION :: tmin, tmax

    CALL fc2d_geoclaw_get_dtopo_interval(tmin, tmax)

    IF (tmin .lt. tmax) THEN
        !! dtopo time is a finite interval : dtopo interval is [tmin, tmax]
        tau = (t - tmin)/(tmax-tmin)
        if (tau .ge. 0 .and. tau .le. 1) then
            fc2d_geoclaw_check_dtopotime = .true.
        else
            fc2d_geoclaw_check_dtopotime = .false.            
        endif
    ELSE
        !! tmin == tmax : dtopo interval is [tmin, \infty]
        IF (t .ge. tmax) THEN
            !! t is in [tmin, \infty]
            tau = 1.d99
            fc2d_geoclaw_check_dtopotime = .true.
        ELSE
            !! We haven't yet reached the dtopo interval
            tau = -1
            fc2d_geoclaw_check_dtopotime = .false.
        ENDIF
    ENDIF

    RETURN

END FUNCTION fc2d_geoclaw_check_dtopotime





