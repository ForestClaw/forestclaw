SUBROUTINE user5_b4step2_manifold(mx,my,mbc, &
     dx,dy,time,maux,aux,blockno,xd,yd,zd)
  IMPLICIT NONE

  INTEGER mbc, mx, my, maux, blockno
  DOUBLE PRECISION dx, dy, time
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  DOUBLE PRECISION xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
  DOUBLE PRECISION yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
  DOUBLE PRECISION zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

  CALL user5_set_velocity_manifold(mx,my,mbc,dx,dy, &
       time,blockno,xd,yd,zd,aux,maux)

END SUBROUTINE user5_b4step2_manifold
