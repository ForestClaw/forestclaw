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

    b = -1.d0
    slope = 0
    d2xzb = 0
END SUBROUTINE bathy_complete
