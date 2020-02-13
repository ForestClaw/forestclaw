      subroutine annulus46_setaux(mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux,blockno,
     &      area, xd,yd,zd, edgelengths,xnormals,ynormals,
     &      xtangents, ytangents, surfnormals)
      implicit none

      integer mbc, mx,my, meqn, maux
      integer blockno
      double precision dx,dy, xlower, ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      integer i,j, k
      double precision dxdy,t

      include "metric_terms.i"

c     # ----------------------------------------------------------------
c     # Conservative form (cell-centered velocities)
c     # 1       Capacity
c     # 2-3     Cell-centered velocities 
c     # 4-11    Normals at all four faces
c     # 12-15   Edge lengths
c     # ----------------------------------------------------------------      

      dxdy = dx*dy

c     # Capacity : entry (1)
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(i,j,1) = area(i,j)/dxdy
         enddo
      enddo


c     # Center velocities : entries (2-5)      
      t = 0
      call annulus46_set_center_velocities(mx,my,mbc,dx,dy,
     &          blockno,xlower,ylower, t,
     &          edgelengths,xnormals,ynormals,
     &          xtangents, ytangents, surfnormals,
     &          aux, maux)

c     # Needed to scale speeds in Riemann solver when using
c     # cell-centered velocities
      do i = 1-mbc,mx+mbc
          do j = 1-mbc,my+mbc
c             # x-face and y-face edge lengths (6,7)      
              aux(i,j,12) = edgelengths(i,j,1)/dy
              aux(i,j,13) = edgelengths(i+1,j,1)/dy
              aux(i,j,14) = edgelengths(i,j,2)/dx
              aux(i,j,15) = edgelengths(i,j+1,2)/dx
          end do
      end do


      return
      end


      subroutine annulus46_set_center_velocities(mx,my,mbc,
     &          dx,dy,blockno,xlower,ylower,t,
     &          edgelengths,xnormals,ynormals,
     &          xtangents, ytangents, surfnormals,
     &          aux, maux)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy, xlower,ylower, t
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision xc,yc,xc1,yc1,zc1, x,y
      double precision nv(3), vel(3)

      double precision nl(3), nr(3), nb(3), nt(3)

      integer*8 cont, get_context

      integer i,j, k
      double precision wl(3), wb(3)

      include "metric_terms.i"

      cont = get_context()

c     # Cell-centered velocities : entries (4,5,6) 
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

c           # This is not the annulus mapping, but rather maps 
c           # the brick to a unit square      
            call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

            call annulus_center_velocity(xc1,yc1,t,vel)

            aux(i,j,2) = vel(1)
            aux(i,j,3) = vel(2)

            do k = 1,3
                nl(k)  = xnormals(i,  j,  k)
                nr(k)  = xnormals(i+1,j,  k)
                nb(k)  = ynormals(i,  j,  k)
                nt(k)  = ynormals(i,  j+1,k)
            enddo

            aux(i,j,4) = nl(1)
            aux(i,j,5) = nl(2)
            aux(i,j,6) = nr(1)
            aux(i,j,7) = nr(2)
            aux(i,j,8) = nb(1)
            aux(i,j,9) = nb(2)
            aux(i,j,10) = nt(1)
            aux(i,j,11) = nt(2)
         enddo
      enddo

      end

      subroutine setaux_cross(u,v,w)
      implicit none

      double precision u(3), v(3), w(3)

      w(1) =  u(2)*v(3) - u(3)*v(2)
      w(2) = -u(1)*v(3) + u(3)*v(1)
      w(3) =  u(1)*v(2) - u(2)*v(1)

      end




