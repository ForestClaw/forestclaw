      subroutine transport5_setaux_manifold(blockno, mx,my,mbc,
     &      xlower,ylower,dx,dy, area, edgelengths, aux, maux)
      implicit none

      integer mbc, mx,my, maux, blockno
      double precision dx,dy, xlower, ylower
      double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision edgelengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j
      double precision dxdy

c     # ----------------------------------------------------------------
c     # Conservative equation (edge velocities)
c     # 1      capacity
c     # 2-3    Edge velocities at left/right x faces
c     # 4-5    Edge velocities at top/bottom y faces
c     # 6-7    Edge lengths (x-faces, y-faces)
c     # 
c     # edge velocities are cell-centered velocities projected onto
c     # edge normals at each of the four faces.
c     # ----------------------------------------------------------------


      dxdy = dx*dy

c     # Capacity : entry (1)
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(1,i,j) = area(i,j)/dxdy
         enddo
      enddo

c     # Needed to scale speeds in Riemann solver when using
c     # cell-centered velocities
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-face and y-face edge lengths (6,7)      
            aux(6,i,j) = edgelengths(i,j,1)/dy
            aux(7,i,j) = edgelengths(i,j,2)/dx
         enddo
      enddo

      return
      end


C       subroutine square5_set_velocities(blockno, mx,my,mbc,
C      &        dx,dy,xlower,ylower, t, xnormals,ynormals,
C      &        surfnormals, aux, maux)
C       implicit none

C       integer mx,my,mbc,maux,blockno
C       double precision dx,dy, xlower,ylower, t
C       double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
C       double precision     xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
C       double precision     ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
C       double precision surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)

C       double precision xc1,yc1,nv(3), vel(3), vdotn, map_dot

C       double precision nl(3), nr(3), nb(3), nt(3)
C       double precision urrot, ulrot, ubrot, utrot

C       integer i,j, k

C c      include "metric_terms.i"

C c     # Cell-centered velocities : entries (4,5,6) 
C       do i = 1-mbc,mx+mbc
C          do j = 1-mbc,my+mbc

C             xc1 = aux(8,i,j)
C             yc1 = aux(9,i,j)

C             call velocity_components_cart(xc1,yc1,t, vel)

C c           # subtract out normal components
C             do k = 1,3
C                 nv(k) = surfnormals(i,j,k)
C             enddo

C             vdotn = map_dot(vel,nv)

C c           # Subtract out component in the normal direction
C             do k = 1,3
C                 vel(k) = vel(k) - vdotn*nv(k)
C             end do

C             do k = 1,3
C                 nl(k)  = xnormals(i,  j,  k)
C                 nr(k)  = xnormals(i+1,j,  k)
C                 nb(k)  = ynormals(i,  j,  k)
C                 nt(k)  = ynormals(i,  j+1,k)
C             enddo

C             ulrot = map_dot(nl,vel)
C             urrot = map_dot(nr,vel)
C             ubrot = map_dot(nb,vel)
C             utrot = map_dot(nt,vel)

C             aux(2,i,j) = ulrot
C             aux(3,i,j) = urrot
C             aux(4,i,j) = ubrot
C             aux(5,i,j) = utrot
C          enddo
C       enddo

C       end



