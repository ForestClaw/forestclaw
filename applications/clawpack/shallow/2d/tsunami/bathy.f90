DOUBLE PRECISION FUNCTION bathy(x,y)
  IMPLICIT NONE

  DOUBLE PRECISION x,y

  DOUBLE PRECISION b, grad(2), d2xzb, d2yzb, d2xyzb

  call sgn_fort_bathy_complete(x,y,b,grad,d2xzb,d2yzb, d2xyzb)

  bathy = b

END FUNCTION bathy

SUBROUTINE sgn_fort_bathy_complete(xc,yc,b,grad,d2xzb,d2yzb,d2xyzb)
    IMPLICIT NONE

    DOUBLE PRECISION xc,yc,b,grad(2),d2xzb, d2yzb, d2xyzb

    b = -1.d0
    grad = 0
    d2xzb = 0
    d2yzb = 0
    d2xyzb = 0
    
END SUBROUTINE sgn_fort_bathy_complete
