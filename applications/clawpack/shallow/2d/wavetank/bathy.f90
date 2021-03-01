DOUBLE PRECISION FUNCTION bathy(x)
  IMPLICIT NONE

  DOUBLE PRECISION x

  double precision b, slope, d2xzb

  call bathy_complete(x,b,slope,d2xzb)

  bathy = b

END FUNCTION bathy

SUBROUTINE bathy_complete(xc,b,slope,d2xzb)
    IMPLICIT NONE

    DOUBLE PRECISION xc,b,slope,d2xzb

    DOUBLE PRECISION xlower, ylower, x
    DOUBLE PRECISION xvals(5), yvals(5)

    DATA xvals /-160.8, -125.5, -90, 0, 40/
    DATA yvals /-4, -4, -0.45, 0, 2/

    xlower = -160.8
    xupper = 40.

!!  # Assume ghost cells use values for first/last end points. 
    if (xc .lt. xlower) then
        x = xlower
    else if (xc .gt. xupper) then
        x = xupper
    endif


    !! From RJL "make_celledges.txt"
    !!  xzpairs = [(-160.8, -4),       # left edge
    !!             (-125.5, -4),       # start of first slope
    !!             (-90,    -0.45),    # start of beach
    !!             (  0,     0),       # shore
    !!             ( 40,     2)]       # right edge

    !! Second derivative is always 0
    d2xzb = 0
    do i = 1,4
        if (xvals(i) .le. x .and. x .le. xvals(i+1)) then
            slope = (yvals(i+1) - yvals(i))/(xvals(i+1) - xvals(i))
            b = yvals(i) + slope*(x - xvals(i))
            return
        endif
    end do

END SUBROUTINE bathy_complete
