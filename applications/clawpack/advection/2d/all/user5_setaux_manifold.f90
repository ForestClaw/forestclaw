SUBROUTINE user5_setaux_manifold(mbc,mx,my, xlower,ylower,dx,dy, &
     maux,aux,blockno,xd,yd,zd,area)
  IMPLICIT NONE

  INTEGER mbc, mx,my, maux
  INTEGER blockno
  DOUBLE PRECISION dx,dy, xlower, ylower
  DOUBLE PRECISION  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  DOUBLE PRECISION xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
  DOUBLE PRECISION yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
  DOUBLE PRECISION zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

  DOUBLE PRECISION area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

  INTEGER i,j
  DOUBLE PRECISION dxdy, t

  dxdy = dx*dy

  DO i = 1-mbc,mx+mbc
     DO j = 1-mbc,my+mbc
        aux(1,i,j) = area(i,j)/dxdy
     ENDDO
  ENDDO

  t = 0
  CALL user5_set_velocity_manifold(mx,my,mbc,dx,dy, &
       t,blockno,xd,yd,zd,aux,maux)

  RETURN
END SUBROUTINE user5_setaux_manifold


SUBROUTINE user5_set_velocity_manifold(mx,my,mbc, &
     dx,dy,t,blockno,xd,yd,zd,aux,maux)
  IMPLICIT NONE

  INTEGER mx,my,mbc,maux, blockno
  DOUBLE PRECISION dx,dy, t

  DOUBLE PRECISION xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
  DOUBLE PRECISION yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
  DOUBLE PRECISION zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

  DOUBLE PRECISION xd1(3),xd2(3)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j
  DOUBLE PRECISION vn
  LOGICAL ispillowsphere

  DO i = 1-mbc,mx+mbc
     DO j = 1-mbc,my+mbc
        !! # x-faces
        xd1(1) = xd(i,j+1)
        xd1(2) = yd(i,j+1)
        xd1(3) = zd(i,j+1)

        xd2(1) = xd(i,j)
        xd2(2) = yd(i,j)
        xd2(3) = zd(i,j)

        CALL get_psi_vel(xd1,xd2,dy,vn,t)
        IF (ispillowsphere()) THEN
           IF (blockno == 1) THEN
              vn = -vn
           ENDIF
        ENDIF
        aux(2,i,j) = vn
     ENDDO
  ENDDO

  DO j = 1-mbc,my+mbc
     DO i = 1-mbc,mx+mbc
        !! # y-faces
        xd1(1) = xd(i+1,j)
        xd1(2) = yd(i+1,j)
        xd1(3) = zd(i+1,j)

        xd2(1) = xd(i,j)
        xd2(2) = yd(i,j)
        xd2(3) = zd(i,j)

        CALL get_psi_vel(xd1,xd2,dx,vn,t)
        IF (ispillowsphere()) THEN
           IF (blockno == 1) THEN
              vn = -vn
           ENDIF
        ENDIF

        aux(3,i,j) = -vn
     ENDDO
  ENDDO

END SUBROUTINE user5_set_velocity_manifold
