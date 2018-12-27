SUBROUTINE torus5_setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
  IMPLICIT NONE

  INTEGER mbc, mx,my, maux
  DOUBLE PRECISION dx,dy, xlower, ylower
  DOUBLE PRECISION  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER blockno, fc2d_clawpack5_get_block
  double precision t

  blockno = fc2d_clawpack5_get_block()

  CALL torus5_velocity_psi_comp(mx,my,mbc,dx,dy, &
       blockno,xlower,ylower,aux,maux)

  RETURN
END SUBROUTINE torus5_setaux


SUBROUTINE torus5_velocity_psi_comp(mx,my,mbc, &
     dx,dy,blockno,xlower,ylower,aux,maux)
  IMPLICIT NONE

  INTEGER mx,my,mbc,maux,blockno
  DOUBLE PRECISION dx,dy, xlower,ylower

  DOUBLE PRECISION xd1(2),xd2(2), t
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j
  DOUBLE PRECISION vn

  t = 0.d0 !! not used

  DO i = 1-mbc,mx+mbc
     DO j = 1-mbc,my+mbc
        !! # x-faces
        xd1(1) = xlower + (i-1)*dx
        xd1(2) = ylower + j*dy

        xd2(1) = xlower + (i-1)*dx
        xd2(2) = ylower + (j-1)*dy

        CALL torus_edge_velocity(blockno,xd1,xd2,dy,vn,t)
        aux(2,i,j) = vn
     ENDDO
  ENDDO

  DO j = 1-mbc,my+mbc
     DO i = 1-mbc,mx+mbc
        !! # y-faces
        xd1(1) = xlower + i*dx
        xd1(2) = ylower + (j-1)*dy

        xd2(1) = xlower + (i-1)*dx
        xd2(2) = ylower + (j-1)*dy

        CALL torus_edge_velocity(blockno,xd1,xd2,dx,vn,t)
        aux(3,i,j) = -vn
     ENDDO
  ENDDO

END SUBROUTINE torus5_velocity_psi_comp
