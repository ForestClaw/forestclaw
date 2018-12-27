c      subroutine torus46_setaux(maxmx, maxmy, mbc,mx,my,
c     &      xlower,ylower, dx,dy, maux,aux)
c      implicit none
c
c      integer maxmx, maxmy, mbc, mx,my, meqn, maux
c      double precision dx,dy, xlower, ylower, z
c      double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)
c
c      integer i,j, blockno, fc2d_clawpack46_get_block
c      double precision dxdy, t, xll, yll, psi
c
c      integer example
c
c      common /excomm_example/ example
c
c
c      blockno = fc2d_clawpack46_get_block()
c
c
c      if (example <= 1) then
c         call torus46_velocity_psi_comp(mx,my,mbc,dx,dy,
c     &       blockno,xlower,ylower,aux,maux)
c      endif
c      
c
c      return
c      end



c      subroutine torus46_set_center_velocities(mx,my,mbc,
c     &      dx,dy,blockno,xlower,ylower,aux,maux)
c      implicit none
c
c      integer mx,my,mbc,maux,blockno
c      double precision dx,dy, xlower,ylower
c
c      double precision xc,yc, t
c      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
c
c      integer i,j
c      double precision vn
c
c      t = 0.d0  !! not used
c      do i = 1-mbc,mx+mbc
c         do j = 1-mbc,my+mbc
cc           # x-faces
c            xc = xlower + (i-0.5)*dx
c            yc = ylower + (j-0.5)*dy
c
c            torus_center_velocity(blockno,xc,yc,t, u, v)
c            aux(i,j,2) = u
c            aux(i,j,3) = v
c         enddo
c      enddo
c
c      end

