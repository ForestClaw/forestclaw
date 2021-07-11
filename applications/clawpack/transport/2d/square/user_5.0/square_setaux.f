      subroutine square5_setaux(blockno, mx,my,mbc,
     &      xlower,ylower,dx,dy, area, edgelengths, 
     &      xp,yp,zp, aux, maux)
      implicit none

      integer mbc, mx,my, maux, blockno
      double precision dx,dy, xlower, ylower
      double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision edgelengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)
      double precision   xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j
      double precision dxdy, xc1, yc1

c     # ----------------------------------------------------------------
c     # Color equation (edge velocities)
c     # 1      capacity
c     # 2-3    Edge velocities at left/right x faces
c     # 4-5    Edge velocities at top/bottom y faces
c     # 6-7    Edge lengths (x-faces, y-faces)
c     3 8-9    Spherical coordinates
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

      do i = 1-mbc,mx+mbc
          do j = 1-mbc,my+mbc
c             # Map to spherical coordinates in [0,1]x[0,1]
              call map2comp(xp(i,j),yp(i,j),zp(i,j),xc1,yc1)

              aux(8,i,j) = xc1
              aux(9,i,j) = yc1
          end do
      end do

      return
      end


      subroutine square5_set_velocities(blockno, mx,my,mbc,
     &        dx,dy,xlower,ylower, t, xnormals,ynormals,
     &        surfnormals, aux, maux)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy, xlower,ylower, t
      double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision     xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision     ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)

      double precision xc1,yc1,nv(3), vel(3), vdotn, map_dot

      double precision nl(3), nr(3), nb(3), nt(3)
      double precision urrot, ulrot, ubrot, utrot

      integer i,j, k

c      include "metric_terms.i"

c     # Cell-centered velocities : entries (4,5,6) 
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc

            xc1 = aux(8,i,j)
            yc1 = aux(9,i,j)

            call velocity_components_cart(xc1,yc1,t, vel)

c           # subtract out normal components
            do k = 1,3
                nv(k) = surfnormals(i,j,k)
            enddo

            vdotn = map_dot(vel,nv)

c           # Subtract out component in the normal direction
            do k = 1,3
                vel(k) = vel(k) - vdotn*nv(k)
            end do

            do k = 1,3
                nl(k)  = xnormals(i,  j,  k)
                nr(k)  = xnormals(i+1,j,  k)
                nb(k)  = ynormals(i,  j,  k)
                nt(k)  = ynormals(i,  j+1,k)
            enddo

            ulrot = map_dot(nl,vel)
            urrot = map_dot(nr,vel)
            ubrot = map_dot(nb,vel)
            utrot = map_dot(nt,vel)

            aux(2,i,j) = ulrot
            aux(3,i,j) = urrot
            aux(4,i,j) = ubrot
            aux(5,i,j) = utrot
         enddo
      enddo

      end



