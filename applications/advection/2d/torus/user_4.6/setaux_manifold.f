      subroutine setaux_manifold(mbc,mx,my,
     &      xlower,ylower,dx,dy,
     &      maux,aux,xp,yp,zp,xd,yd,zd,
     &      xnormals,ynormals,edge_lengths,area)
      implicit none

      integer mbc, mx,my, meqn, maux
      double precision dx,dy, xlower, ylower, z
      double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      integer i,j
      double precision dxdy, t, xll, yll, psi

      include 'metric_terms.i'

      dxdy = dx*dy
      t = 0
      z = 0

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(i,j,1) = area(i,j)/dxdy
         enddo
      enddo

      call compute_velocity_psi_comp(mx,my,mbc,dx,dy,
     &      t,xlower,ylower,aux,maux)

c     # Compute the velocity in Cartesian coordinates;
c     # project onto edge normals and scale. This
c     # may not be strictly speaking conservative.

c      call compute_velocity_ec(mx,my,mbc,meqn,
c     &      xlower,ylower,dx,dy,t,xnormals,ynormals,
c     &      edge_lengths,aux,maux)


      return
      end


      subroutine compute_velocity_ec(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy,t,xnormals,ynormals,
     &      edge_lengths,aux,maux)
      implicit none

      integer mx,my,mbc,maux, meqn
      double precision xlower,ylower,dx,dy, t

      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j,m, mq
      double precision vel(3), xc,yc,xe,ye,ze,v,nv(3)
      double precision ds

      include 'metric_terms.i'

c     # Get edge centered velocities.

      do i = 2-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # Do x-edges
            xe = xlower + (i-1)*dx
            yc = ylower + (j-0.5)*dy

            call get_vel(xe,yc,vel,t)

            v = 0
            do m = 1,3
               nv(m) = xnormals(i,j,m)
               v = v + nv(m)*vel(m)
            enddo
c           # Solution is already scaled by lengths
            ds = edge_lengths(i,j,1)
            aux(i,j,2) = ds*v/dy
         enddo
      enddo

      do j = 2-mbc,my+mbc
         do i = 1-mbc,mx+mbc
c           # Do y-face
            xc = xlower + (i-0.5)*dx
            ye = ylower + (j-1)*dy

            call get_vel(xc,ye,vel,t)

            v = 0
            do m = 1,3
               nv(m) = ynormals(i,j,m)
               v = v + nv(m)*vel(m)
            enddo
            ds = edge_lengths(i,j,2)
            aux(i,j,3) = ds*v/dx
         enddo
      enddo

      end


      subroutine compute_velocity_psi_comp(mx,my,mbc,
     &      dx,dy,t,xlower,ylower,aux,maux)
      implicit none

      integer mx,my,mbc,maux
      double precision dx,dy, t, xlower,ylower

      double precision xd1(2),xd2(2)
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j
      double precision vn

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-faces
            xd1(1) = xlower + (i-1)*dx
            xd1(2) = ylower + j*dy

            xd2(1) = xlower + (i-1)*dx
            xd2(2) = ylower + (j-1)*dy

            call get_vel_psi_comp(xd1,xd2,dy,vn,t)
            aux(i,j,2) = vn
         enddo
      enddo

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
c           # y-faces
            xd1(1) = xlower + i*dx
            xd1(2) = ylower + (j-1)*dy

            xd2(1) = xlower + (i-1)*dx
            xd2(2) = ylower + (j-1)*dy

            call get_vel_psi_comp(xd1,xd2,dx,vn,t)
            aux(i,j,3) = -vn
         enddo
      enddo

      end

      subroutine get_vel_psi_comp(xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(2),xd2(2), ds, vn, psi,t

      vn = (psi(xd1(1),xd1(2),t) -
     &      psi(xd2(1),xd2(2),t))/ds

      end
