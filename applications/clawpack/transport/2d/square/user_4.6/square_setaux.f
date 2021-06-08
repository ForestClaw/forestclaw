      subroutine square_setaux(blockno, mx,my,mbc,
     &      xlower,ylower,dx,dy, area, edgelengths,
     &      xnormals,ynormals, surfnormals, aux, maux)
      implicit none

      integer mbc, mx,my, meqn, maux
      integer blockno
      double precision dx,dy, xlower, ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      integer i,j, k
      double precision dxdy

      include "fclaw2d_metric_terms.i"

c     # ----------------------------------------------------------------
c     # Color equation (edge velocities)
c     # 1      capacity
c     # 2-3    Edge velocities at left/right x faces
c     # 4-5    Edge velocities at top/bottom y faces
c     # 6-7    Edge lengths (x-faces, y-faces)
c     # ----------------------------------------------------------------


      dxdy = dx*dy

c     # Capacity : entry (1)
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(i,j,1) = area(i,j)/dxdy
         enddo
      enddo

c     # Center velocities : entries (2-5)      
      call square_set_center_velocities(mx,my,mbc,dx,dy,
     &          blockno,xlower,ylower,
     &          edgelengths,xnormals,ynormals,surfnormals,
     &          aux, maux)

c     # Needed to scale speeds in Riemann solver when using
c     # cell-centered velocities
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-face and y-face edge lengths (6,7)      
            aux(i,j,6) = edgelengths(i,j,1)/dy
            aux(i,j,7) = edgelengths(i,j,2)/dx
         enddo
      enddo

      return
      end


      subroutine square_set_center_velocities(mx,my,mbc,
     &          dx,dy,blockno,xlower,ylower,
     &          edgelengths,xnormals,ynormals,surfnormals,
     &          aux, maux)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy, xlower,ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision xc,yc
      double precision xc1,yc1,zc1, nv(3), vel(3), vdotn, map_dot

      double precision nl(3), nr(3), nb(3), nt(3)
      double precision urrot, ulrot, ubrot, utrot
      double precision x,y

      integer*8 cont, get_context

      integer i,j, k

      include "fclaw2d_metric_terms.i"

      cont = get_context()

c     # Cell-centered velocities : entries (4,5,6) 
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

c           # This is not the square mapping, but rather maps the brick to
c           # a unit square      
            call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

            call square_center_velocity(blockno, xc1,yc1,vel)

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

c            ulrot = nl(1)*vel(1) + nl(2)*vel(2) + nl(3)*vel(3)
c            urrot = nr(1)*vel(1) + nr(2)*vel(2) + nr(3)*vel(3)
c            ubrot = nb(1)*vel(1) + nb(2)*vel(2) + nb(3)*vel(3)
c            utrot = nt(1)*vel(1) + nt(2)*vel(2) + nt(3)*vel(3)

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



