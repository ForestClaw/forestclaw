      subroutine torus46_setaux(maxmx, maxmy, mbc,mx,my,
     &      xlower,ylower, dx,dy, maux,aux)
      implicit none

      integer maxmx, maxmy, mbc, mx,my, meqn, maux
      double precision dx,dy, xlower, ylower, z
      double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      integer i,j, blockno, fc2d_clawpack46_get_block
      double precision dxdy, t, xll, yll, psi

      blockno = fc2d_clawpack46_get_block()

      call torus46_velocity_psi_comp(mx,my,mbc,dx,dy,
     &      blockno,xlower,ylower,aux,maux)

      return
      end

      subroutine torus46_velocity_psi_comp(mx,my,mbc,
     &      dx,dy,blockno,xlower,ylower,aux,maux)
      implicit none

      integer mx,my,mbc,maux,blockno
      double precision dx,dy, xlower,ylower

      double precision xd1(2),xd2(2), t
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j
      double precision vn

      t = 0.d0  !! not used
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-faces
            xd1(1) = xlower + (i-1)*dx
            xd1(2) = ylower + j*dy

            xd2(1) = xlower + (i-1)*dx
            xd2(2) = ylower + (j-1)*dy

            call get_vel_psi_comp(blockno,xd1,xd2,dy,vn,t)
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

            call get_vel_psi_comp(blockno,xd1,xd2,dx,vn,t)
            aux(i,j,3) = -vn
         enddo
      enddo

      end
