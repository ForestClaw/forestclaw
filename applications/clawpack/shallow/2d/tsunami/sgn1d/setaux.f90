SUBROUTINE setaux(mbc,mx,xlower,dx,maux,aux)

  IMPLICIT NONE
  INTEGER, INTENT(in) :: mbc,mx,maux
  REAL(KIND=8), INTENT(in) :: xlower,dx
  REAL(KIND=8), INTENT(out) ::  aux(maux,1-mbc:mx+mbc)

  INTEGER i, ibc
  DOUBLE PRECISION xc, yc, bathy, b, grad(2), d2xzb, d2yzb, d2xyzb

  if (maux .lt. 3) then
      write(6,*) 'setaux.f : maux must be at least 3'
      stop
  endif

  yc = 0
  DO i = 1-mbc,mx+mbc
     xc = xlower + (i-0.5)*dx
     call sgn_fort_bathy_complete(xc,yc,b,grad,d2xzb,d2yzb,d2xyzb)
     aux(1,i) = b
     !! aux(2,i) = slope
     !! aux(3,i) = d2xzb
  end do

  do i = 2-mbc,mx+1
     aux(2,i) = (aux(1,i+1) - aux(1,i-1))/(2.d0*dx)
     aux(3,i) = (aux(1,i+1) - 2.d0*aux(1,i) + aux(1,i-1))/(dx*dx)
  END DO


END SUBROUTINE SETAUX