      subroutine setaux_manifold(mbc,mx,my,
     &      xlower,ylower,dx,dy, maux,aux,xp,yp,zp,xd,yd,zd,
     &      area,ismanifold)
      implicit none

      integer mbc, mx,my, meqn, maux
      logical ismanifold
      double precision dx,dy, xlower, ylower
      double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, maux)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j
      double precision dxdy, t, sum, pi, exact_area
      logical iscart, issphere, isdisk,debug

      logical fclaw2d_map_is_used
      integer**8 cont
      integer blockno

      common /compi/ pi

c     # keep debug false since only one mpirank should print output
      dxdy = dx*dy

      cont = get_context()

      if (fclaw2d_map_is_used(cont)) then
         sum = 0
         do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc
               aux(i,j,1) = area(i,j)/dxdy
               if ((i .ge. 1 .and. i .le. mx) .and.
     &               (j .ge. 1 .and. j .le. my)) then
                  sum = sum + area(i,j)
               endif
            enddo
         enddo
      endif

      if (debug) then
c        # We don't get an exact area, because of the non-smoothness
c        # across the diagonals
         write(6,100) 'Surface area', sum
         if (isdisk()) then
            exact_area = pi
         elseif (iscart()) then
            exact_area = 4.d0
         elseif (issphere()) then
c           # Hemisphere
            exact_area = 2*pi
         endif
         write(6,100) 'Error', abs(sum - exact_area)

         stop
  100    format(A15,E24.16)
=======
         t = 0
         call compute_velocity_psi(mx,my,mbc,dx,dy,
     &         t,xd,yd,zd,aux,maux)
      else
         call compute_velocity_psi_nomanifold(mx,my,mbc,dx,dy,
     &         xlower,ylower,t,aux,maux)
>>>>>>> Everything should work without computing area,xp,yp,zp,...
      endif

      return
      end


      subroutine compute_velocity_psi(mx,my,mbc,
     &      dx,dy,t,xd,yd,zd,aux,maux)
      implicit none

      integer maxmx, maxmy, mx,my,mbc,maux
      double precision dx,dy, t

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      double precision xd1(3),xd2(3)
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j
      double precision vn

      logical ispillowsphere

      integer blockno, get_block
      blockno = get_block()


      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # x-faces
            xd1(1) = xd(i,j+1)
            xd1(2) = yd(i,j+1)
            xd1(3) = zd(i,j+1)

            xd2(1) = xd(i,j)
            xd2(2) = yd(i,j)
            xd2(3) = zd(i,j)

            call get_vel_psi(xd1,xd2,dy,vn,t)
            if (ispillowsphere()) then
               if (blockno .eq. 1) then
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

            call get_vel_psi(xd1,xd2,dx,vn,t)
            if (ispillowsphere()) then
               if (blockno .eq. 1) then
                  vn = -vn
               endif
            endif

            aux(i,j,3) = -vn
         enddo
      enddo

      end


      subroutine compute_velocity_psi_nomanifold(mx,my,mbc,
     &      xlower,ylower,dx,dy,t,aux,maux)
      implicit none

      integer mx,my,mbc,maux
      double precision dx,dy, t, xlower,ylower


      double precision xd1(3),xd2(3)
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j
      double precision vn

      integer blockno, get_block
      blockno = get_block()


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

      double precision xd1(3),xd2(3), ds, vn, psi,t

      vn = (psi(xd1(1),xd1(2),xd1(3)) -
     &      psi(xd2(1),xd2(2),xd2(3)))/ds

      end
