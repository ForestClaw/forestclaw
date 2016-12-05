      subroutine user46_setaux_manifold(mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux,blockno,
     &      xd,yd,zd,area)
      implicit none

      integer mbc, mx,my, meqn, maux
      integer blockno
      double precision dx,dy, xlower, ylower
      double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j
      double precision dxdy, t, x

      dxdy = dx*dy

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            aux(i,j,1) = area(i,j)/dxdy
         enddo
      enddo

      t = 0
      call user46_set_velocity_manifold(mx,my,mbc,
     &      dx,dy,t,blockno,xd,yd,zd,aux,maux)

      return
      end

      subroutine user46_set_velocity_manifold(mx,my,mbc,
     &      dx,dy,t,blockno,xd,yd,zd,aux,maux)
      implicit none

      integer mx,my,mbc,maux, blockno
      double precision dx,dy, t

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      double precision xd1(3),xd2(3)
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j
      double precision vn
      logical ispillowsphere

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-faces
            xd1(1) = xd(i,j+1)
            xd1(2) = yd(i,j+1)
            xd1(3) = zd(i,j+1)

            xd2(1) = xd(i,j)
            xd2(2) = yd(i,j)
            xd2(3) = zd(i,j)

            call get_psi_vel(xd1,xd2,dy,vn,t)
            if (ispillowsphere()) then
               if (blockno == 1) then
                  vn = -vn
               endif
            endif
            aux(i,j,2) = vn
         enddo
      enddo

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
c           # y-faces
            xd1(1) = xd(i+1,j)
            xd1(2) = yd(i+1,j)
            xd1(3) = zd(i+1,j)

            xd2(1) = xd(i,j)
            xd2(2) = yd(i,j)
            xd2(3) = zd(i,j)

            call get_psi_vel(xd1,xd2,dx,vn,t)
            if (ispillowsphere()) then
               if (blockno == 1) then
                  vn = -vn
               endif
            endif

            aux(i,j,3) = -vn
         enddo
      enddo
      end
