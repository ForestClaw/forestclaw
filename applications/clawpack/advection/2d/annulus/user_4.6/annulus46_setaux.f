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
      double precision dxdy

      integer example
      common /example_comm/ example  

      integer color_equation
      common /eqn_comm/ color_equation      

      include "metric_terms.i"

c     # ----------------------------------------------------------------
c     # Color equation (edge velocities)
c     # 1      capacity
c     # 2-3    Edge velocities
c     #
c     # Conservative form (cell-centered velocities)
c     # 2-5    Cell-centered velocities projected onto four edge normals
c     # 6-7    Edge lengths (x-face, y-face)
c     # ----------------------------------------------------------------

      dxdy = dx*dy

c     # Capacity : entry (1)
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(i,j,1) = area(i,j)/dxdy
         enddo
      enddo

      if (color_equation .eq. 1) then
c         # Edge velocities using a streamfunction : entries (2-3)      
          call annulus46_set_edge_velocities(mx,my,mbc,dx,dy,
     &          blockno,xd,yd,aux,maux)
      else
c         # Center velocities : entries (2-5)      
          call annulus46_set_center_velocities(mx,my,mbc,dx,dy,
     &          blockno,xlower,ylower,
     &          edgelengths,xnormals,ynormals,
     &          xtangents, ytangents, surfnormals,
     &          aux, maux)

c         # Needed to scale speeds in Riemann solver when using
c         # cell-centered velocities
          do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc
c               # x-face and y-face edge lengths (6,7)      
                aux(i,j,6) = edgelengths(i,j,1)/dy
                aux(i,j,7) = edgelengths(i,j,2)/dx
            enddo
          enddo
      endif


      return
      end


      subroutine annulus46_set_edge_velocities(mx,my,mbc,
     &      dx,dy, blockno,xd,yd,aux,maux)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      double precision xd1, yd1, xd2, yd2

      integer i,j
      double precision vn

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc

            xd1 = xd(i,j+1)
            yd1 = yd(i,j+1)
            xd2 = xd(i,j)
            yd2 = yd(i,j)

            call annulus_edge_velocity(xd1,yd1,xd2,yd2,dy,vn)
            aux(i,j,2) = -vn
         enddo
      enddo

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc

            xd1 = xd(i+1,j)
            yd1 = yd(i+1,j)
            xd2 = xd(i,j)
            yd2 = yd(i,j)

            call annulus_edge_velocity(xd1,yd1,xd2,yd2,dx,vn)
            aux(i,j,3) = vn
         enddo
      enddo

      end

      subroutine annulus46_set_center_velocities(mx,my,mbc,
     &          dx,dy,blockno,xlower,ylower,
     &          edgelengths,xnormals,ynormals,
     &          xtangents, ytangents, surfnormals,
     &          aux, maux)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy, xlower,ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision xc,yc,xc1,yc1,zc1, x,y
      double precision nv(3), vel(3), vdotn, annulus_dot

      double precision nl(3), nr(3), nb(3), nt(3)
      double precision tl(3), tr(3), tb(3), tt(3)
      double precision urrot, ulrot, ubrot, utrot

      integer*8 cont, get_context

      integer i,j, k
      double precision wl(3), wb(3)

      integer mapping
      common /mapping_comm/ mapping

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
            call annulus_transform_coordinates(xc1,yc1,x,y,mapping)

            call annulus_center_velocity(x,y,vel)

c           # Subtract out component in the normal direction
            do k = 1,3
                nv(k) = surfnormals(i,j,k)
            enddo

            vdotn = annulus_dot(vel,nv)

            do k = 1,3
                vel(k) = vel(k) - vdotn*nv(k)
            end do

            do k = 1,3
                nl(k)  = xnormals(i,  j,  k)
                nr(k)  = xnormals(i+1,j,  k)
                nb(k)  = ynormals(i,  j,  k)
                nt(k)  = ynormals(i,  j+1,k)
            enddo

            do k = 1,3
                tl(k)  = xtangents(i,  j,  k)
                tr(k)  = xtangents(i+1,j,  k)
                tb(k)  = ytangents(i,  j,  k)
                tt(k)  = ytangents(i,  j+1,k)
            enddo

            ulrot = nl(1)*vel(1) + nl(2)*vel(2) + nl(3)*vel(3)
            urrot = nr(1)*vel(1) + nr(2)*vel(2) + nr(3)*vel(3)
            ubrot = nb(1)*vel(1) + nb(2)*vel(2) + nb(3)*vel(3)
            utrot = nt(1)*vel(1) + nt(2)*vel(2) + nt(3)*vel(3)

c            write(6,100) nb(1), nb(2), ubrot
100         format(3F16.8)            

            aux(i,j,2) = ulrot
            aux(i,j,3) = urrot
            aux(i,j,4) = ubrot
            aux(i,j,5) = utrot
         enddo
      enddo
c      stop

      end

      subroutine setaux_cross(u,v,w)
      implicit none

      double precision u(3), v(3), w(3)

      w(1) =  u(2)*v(3) - u(3)*v(2)
      w(2) = -u(1)*v(3) + u(3)*v(1)
      w(3) =  u(1)*v(2) - u(2)*v(1)

      end




