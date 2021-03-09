DOUBLE PRECISION FUNCTION bathy(x,y)
  IMPLICIT NONE

  DOUBLE PRECISION x,y

  double precision b, grad(2), d2xzb, d2yzb, d2xyzb

  call sgn_fort_bathy_complete(x,y,b,grad,d2xzb,d2yzb,d2xyzb)

  bathy = b

END FUNCTION bathy


SUBROUTINE sgn_fort_bathy_complete(xc,yc, b,grad,d2xzb, d2yzb, d2xyzb)
    IMPLICIT NONE

    DOUBLE PRECISION xc,yc, b,grad(2),d2xzb, d2yzb, d2xyzb

    DOUBLE PRECISION xlower, xupper, x
    DOUBLE PRECISION xvals(5), bvals(5)

    DOUBLE PRECISION  grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level


    integer i
    double precision bflat

    DATA xvals /-160.8, -125.5, -90,    0, 40/
    DATA bvals /-4,     -4,      -0.45, 0,  2/

    xlower = -160.8
    xupper = 40.

!!  # Assume ghost cells use values for first/last end points. 
    if (xc .lt. xlower) then
        x = xlower
    else if (xc .gt. xupper) then
        x = xupper
    else
        x = xc
    endif

    !! From RJL "make_celledges.txt"
    !!  xzpairs = [(-160.8, -4),       # left edge
    !!             (-125.5, -4),       # start of first slope
    !!             (-90,    -0.45),    # start of beach
    !!             (  0,     0),       # shore
    !!             ( 40,     2)]       # right edge

    bflat = -4

    !! Second derivative is always 0
    d2xzb = 0
    d2yzb = 0
    d2xyzb = 0
    do i = 1,4
        if (xvals(i) .le. x .and. x .le. xvals(i+1)) then
            grad(1) = 0
            grad(2) = 0
            b = bflat !! bvals(i) + grad(1)*(x - xvals(i))
            return
        endif
    end do

    write(6,*) 'sgn_fort_bathy_complete : No valid x value'
    stop

END SUBROUTINE sgn_fort_bathy_complete
