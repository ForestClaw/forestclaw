c     # This is the "classic" call to setaux.
      subroutine clawpack46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit none

      integer mbc, mx,my, meqn, maux, maxmx, maxmy
      double precision dx,dy, xlower, ylower
      double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      double precision t

      t = 0
      call compute_velocity_psi_nomap(mx,my,mbc,dx,dy,
     &      xlower,ylower,t,aux,maux)

      end

      subroutine compute_velocity_psi_nomap(mx,my,mbc,
     &      dx,dy,xlower,ylower,t,aux,maux)
      implicit none

      integer mx,my,mbc,maux
      double precision dx,dy, t, xlower,ylower


      double precision xd1(3),xd2(3)
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j
      double precision vn

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-faces
            xd1(1) = xlower + (i-1)*dx
            xd1(2) = ylower + (j)*dy
            xd1(3) = 0

            xd2(1) = xlower + (i-1)*dx
            xd2(2) = ylower + (j-1)*dy
            xd2(3) = 0

            call get_vel_psi(xd1,xd2,dy,vn,t)
            aux(i,j,2) = vn
         enddo
      enddo

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
c           # y-faces
            xd1(1) = xlower + i*dx
            xd1(2) = ylower + (j-1)*dy
            xd1(3) = 0

            xd2(1) = xlower + (i-1)*dx
            xd2(2) = ylower + (j-1)*dy
            xd2(3) = 0

            call get_vel_psi(xd1,xd2,dx,vn,t)

            aux(i,j,3) = -vn
         enddo
      enddo

      end

      subroutine get_vel_psi(xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(2),xd2(2), ds, vn, psi_nomap,t
      integer blockno

      vn = (psi_nomap(xd1(1),xd1(2)) -
     &      psi_nomap(xd2(1),xd2(2)))/ds

      end

      double precision function psi_nomap(x,y)
      implicit none

      double precision x,y, revs_per_s
      double precision u0, v0

      double precision u0_comm,v0_comm,revs_comm
      common /comm_velocity/ u0_comm,v0_comm,revs_comm

      u0 = revs_comm*u0_comm
      v0 = revs_comm*v0_comm

      psi_nomap = (-v0*x + u0*y)

      end
