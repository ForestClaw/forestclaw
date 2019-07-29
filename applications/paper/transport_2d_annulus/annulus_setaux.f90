SUBROUTINE setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
  IMPLICIT NONE

  INTEGER mbc, mx,my, meqn, maux
  DOUBLE PRECISION dx,dy, xlower, ylower
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER :: blockno, ghost_only, quadsize
  INTEGER :: i,j, k, m
  DOUBLE PRECISION :: dxdy, xc, yc, nv(3), t, u, v, w
  DOUBLE PRECISION, allocatable, DIMENSION(:) :: quadstore
  DOUBLE PRECISION :: sum, dx0
  INTEGER :: level
  logical :: area_check

  INTEGER color_equation
  COMMON /eqn_comm/ color_equation


  INTEGER maxlevel, rfactor, grid_mx, mi, mj
  COMMON /amr_comm/ maxlevel, rfactor, grid_mx, mi, mj

  INCLUDE "metric_terms.i"

  area_check = .false.

  !! Figure out level of current grid so we can compute areas
  !! in a conservative manner
  dx0 = 1.d0/(mi*grid_mx)
  level = log((dx0+1d-8)/dx)/log(real(rfactor)) + 1

  !!     # ----------------------------------------------------------------
  !!     # Color equation (edge velocities)
  !!     # 1      capacity
  !!     # 2-3    Edge velocities
  !!     #
  !!     # Conservative form (cell-centered velocities)
  !!     # 2-5    Cell-centered velocities projected onto four edge normals
  !!     # 6-7    Edge lengths (x-face, y-face)
  !!     # ----------------------------------------------------------------

  blockno = 0
  ghost_only = 0

  CALL fclaw2d_fort_compute_mesh(mx,my,mbc,&
       xlower,ylower,dx,dy,blockno,xp,yp,zp,xd,yd,zd)

  !! from ForestClaw
  m = rfactor**(maxlevel-level)
  ALLOCATE(quadstore(3*(m+1)*(m+1)))
  quadsize = m
  CALL fclaw2d_fort_compute_area(mx,my,mbc,dx,dy, &
       xlower, ylower, blockno,area, quadsize, quadstore, &
       ghost_only)
  DEALLOCATE(quadstore)

  CALL fclaw2d_fort_compute_normals(mx,my,mbc, &
       xp,yp,zp,xd,yd,zd,xnormals,ynormals)

  CALL fclaw2d_fort_compute_tangents(mx,my,mbc, &
       xd,yd,zd, xtangents,ytangents,edgelengths)

  CALL fclaw2d_fort_compute_surf_normals(mx,my,mbc, &
       xnormals,ynormals,edgelengths,curvature, &
       surfnormals,area)

  t = 0
  dxdy = dx*dy

  !!     # Capacity : entry (1)
  DO i = 1-mbc,mx+mbc
     DO j = 1-mbc,my+mbc
        aux(1,i,j) = area(i,j)/dxdy
     END DO
  END DO

  if (area_check) then
    sum = 0
    do i = 1,mx
        do j = 1,my
            sum = sum + area(i,j)
        end do
    end do
    write(6,101) level,sum
101 format(I5,F24.16) 
endif

  IF (color_equation .EQ. 1) THEN
     !! # Edge velocities using a streamfunction : entries (2-3)
     CALL annulus46_set_edge_velocities(mx,my,mbc,dx,dy, &
          blockno,xd,yd,aux,maux)
  ELSE
     !! # Center velocities : entries (2-5)
     CALL annulus46_set_center_velocities(mx,my,mbc,dx,dy, &
          blockno,xlower,ylower, &
          edgelengths,xnormals,ynormals,surfnormals, &
          aux, maux)

     !! # Needed to scale speeds in Riemann solver when using
     !! # cell-centered velocities
     DO i = 1-mbc,mx+mbc
        DO j = 1-mbc,my+mbc
           !! Need to store all four edge lengths in each cell, since they are all
           !! needed in QAD fix.            
           aux(6,i,j) = edgelengths(i,  j,  1)/dy
           aux(7,i,j) = edgelengths(i+1,j,  1)/dy
           aux(8,i,j) = edgelengths(i,  j,  2)/dx
           aux(9,i,j) = edgelengths(i,  j+1,2)/dx
        ENDDO
     ENDDO
  ENDIF


END SUBROUTINE setaux


SUBROUTINE annulus46_set_edge_velocities(mx,my,mbc, &
     dx,dy, blockno,xd,yd,aux,maux)
  IMPLICIT NONE

  INTEGER mx,my,mbc,maux,blockno
  DOUBLE PRECISION dx,dy
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  DOUBLE PRECISION xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
  DOUBLE PRECISION yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

  DOUBLE PRECISION xd1, yd1, xd2, yd2

  INTEGER i,j
  DOUBLE PRECISION vn

  DO i = 1-mbc,mx+mbc
     DO j = 1-mbc,my+mbc

        xd1 = xd(i,j+1)
        yd1 = yd(i,j+1)
        xd2 = xd(i,j)
        yd2 = yd(i,j)

        CALL annulus_edge_velocity(xd1,yd1,xd2,yd2,dy,vn)
        aux(2,i,j) = -vn
     ENDDO
  ENDDO

  DO j = 1-mbc,my+mbc
     DO i = 1-mbc,mx+mbc

        xd1 = xd(i+1,j)
        yd1 = yd(i+1,j)

        xd2 = xd(i,j)
        yd2 = yd(i,j)

        CALL annulus_edge_velocity(xd1,yd1,xd2,yd2,dx,vn)
        aux(3,i,j) = vn
     ENDDO
  ENDDO

END SUBROUTINE annulus46_set_edge_velocities


SUBROUTINE annulus46_set_center_velocities(mx,my,mbc, &
     dx,dy,blockno,xlower,ylower, &
     edgelengths,xnormals,ynormals,surfnormals, &
     aux, maux)
  IMPLICIT NONE

  INTEGER mx,my,mbc,maux,blockno
  DOUBLE PRECISION dx,dy, xlower,ylower
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  DOUBLE PRECISION xc,yc,xc1,yc1,zc1, x,y
  DOUBLE PRECISION nv(3), vel(3), vdotn, annulus_dot

  DOUBLE PRECISION nl(3), nr(3), nb(3), nt(3)
  DOUBLE PRECISION urrot, ulrot, ubrot, utrot

  INTEGER i,j, k

  DOUBLE PRECISION beta
  COMMON /annulus_comm/ beta

  INTEGER mapping
  COMMON /mapping_comm/ mapping

  INCLUDE "metric_terms.i"


  !! # Cell-centered velocities : entries (4,5,6)
  DO i = 1-mbc,mx+mbc
     DO j = 1-mbc,my+mbc
        xc = xlower + (i-0.5)*dx
        yc = ylower + (j-0.5)*dy

        CALL annulus_transform_coordinates(xc,yc,x,y,mapping)
        CALL annulus_center_velocity(x,y,vel)

        !! # Subtract out component in the normal direction
!!        DO k = 1,3
!!           nv(k) = surfnormals(i,j,k)
!!        ENDDO
!!
!!        vdotn = annulus_dot(vel,nv)
!!
!!        DO k = 1,3
!!           vel(k) = vel(k) - vdotn*nv(k)
!!        END DO

        DO k = 1,3
           nl(k)  = xnormals(i,  j,  k)
           nr(k)  = xnormals(i+1,j,  k)
           nb(k)  = ynormals(i,  j,  k)
           nt(k)  = ynormals(i,  j+1,k)
        ENDDO

        ulrot = nl(1)*vel(1) + nl(2)*vel(2) + nl(3)*vel(3)
        urrot = nr(1)*vel(1) + nr(2)*vel(2) + nr(3)*vel(3)
        ubrot = nb(1)*vel(1) + nb(2)*vel(2) + nb(3)*vel(3)
        utrot = nt(1)*vel(1) + nt(2)*vel(2) + nt(3)*vel(3)

        aux(2,i,j) = ulrot
        aux(3,i,j) = urrot
        aux(4,i,j) = ubrot
        aux(5,i,j) = utrot
     ENDDO
  ENDDO
END SUBROUTINE annulus46_set_center_velocities
