      subroutine torus46_setaux(mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux,blockno,
     &      area, edgelengths,xnormals,ynormals,
     &      surfnormals)
      implicit none

      integer mbc, mx,my, meqn, maux
      integer blockno
      double precision dx,dy, xlower, ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      integer i,j, k
      double precision dxdy

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
c            if (i .ge. 1 .and. i .le. mx .and. 
c     &          j .ge. 1 .and. j .le. my) then
c                  write(6,102) i,j,area(i,j)
c            endif
c102         format(2I5,F24.16)            
         enddo
      enddo

c     # Center velocities : entries (2-5)      
      call torus46_set_center_velocities(mx,my,mbc,dx,dy,
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


      subroutine torus46_set_center_velocities(mx,my,mbc,
     &          dx,dy,blockno,xlower,ylower,
     &          edgelengths,xnormals,ynormals,surfnormals,
     &          aux, maux)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy, xlower,ylower
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision xc,yc
      double precision xc1,yc1,zc1, nv(3), vel(3), vdotn, torus_dot

      double precision nl(3), nr(3), nb(3), nt(3)
      double precision urrot, ulrot, ubrot, utrot
      double precision x,y,z

      integer mapping
      common /mapping_comm/ mapping

      integer*8 cont, get_context

      integer i,j, k

      include "metric_terms.i"

      cont = get_context()

c     # Cell-centered velocities : entries (4,5,6) 
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

c           # This is not the torus mapping, but rather maps the brick to
c           # a unit square      
            call fclaw2d_map_brick2c(cont,blockno,xc,yc,x,y,z)

            call torus_center_velocity(x,y,vel)

c           # subtract out normal components
            do k = 1,3
                nv(k) = surfnormals(i,j,k)
            enddo

            vdotn = torus_dot(vel,nv)

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

            ulrot = nl(1)*vel(1) + nl(2)*vel(2) + nl(3)*vel(3)
            urrot = nr(1)*vel(1) + nr(2)*vel(2) + nr(3)*vel(3)
            ubrot = nb(1)*vel(1) + nb(2)*vel(2) + nb(3)*vel(3)
            utrot = nt(1)*vel(1) + nt(2)*vel(2) + nt(3)*vel(3)

            aux(i,j,2) = ulrot
            aux(i,j,3) = urrot
            aux(i,j,4) = ubrot
            aux(i,j,5) = utrot
         enddo
      enddo

      end



