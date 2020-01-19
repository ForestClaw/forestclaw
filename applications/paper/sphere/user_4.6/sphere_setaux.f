      subroutine sphere_setaux(blockno, mx,my,mbc,
     &      xlower,ylower,dx,dy, area, edgelengths, 
     &      xp,yp,zp, aux, maux)
      implicit none

      integer mbc, mx,my, meqn, maux
      integer blockno
      double precision dx,dy, xlower, ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      integer i,j, k
      double precision dxdy, xc1, yc1

      include "metric_terms.i"

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
            aux(i,j,1) = area(i,j)/dxdy
         enddo
      enddo

c     # Needed to scale speeds in Riemann solver when using
c     # cell-centered velocities
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-face and y-face edge lengths (6,7)      
            aux(i,j,6) = edgelengths(i,j,1)/dy
            aux(i,j,7) = edgelengths(i,j,2)/dx
         enddo
      enddo

      do i = 1-mbc,mx+mbc
          do j = 1-mbc,my+mbc
c             # Map to spherical coordinates in [0,1]x[0,1]
              call map2comp(xp(i,j),yp(i,j),zp(i,j),xc1,yc1)

              aux(i,j,8) = xc1
              aux(i,j,9) = yc1
          end do
      end do

      return
      end


      subroutine sphere_set_velocities(blockno, mx,my,mbc,
     &        dx,dy,xlower,ylower, t, xnormals,ynormals,
     &        surfnormals, aux, maux)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy, xlower,ylower, t
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision xc,yc
      double precision xc1,yc1,zc1, nv(3), vel(3), vdotn, map_dot

      double precision nl(3), nr(3), nb(3), nt(3)
      double precision urrot, ulrot, ubrot, utrot

      integer curvature_correction
      common /conservation_com/ curvature_correction

c      integer*8 cont, get_context

      integer i,j, k

      include "metric_terms.i"

c      cont = get_context()

c     # Cell-centered velocities : entries (4,5,6) 
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc

            xc1 = aux(i,j,8)
            yc1 = aux(i,j,9)

            call velocity_components_cart(xc1,yc1,t, vel)

c           # subtract out normal components
            do k = 1,3
                nv(k) = surfnormals(i,j,k)
            enddo

c           # Subtract out component in the normal direction
            if (curvature_correction .ne. 0) then
                vdotn = map_dot(vel,nv)
                do k = 1,3
                    vel(k) = vel(k) - vdotn*nv(k)
                end do
            endif

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

            aux(i,j,2) = ulrot
            aux(i,j,3) = urrot
            aux(i,j,4) = ubrot
            aux(i,j,5) = utrot
         enddo
      enddo

      end



