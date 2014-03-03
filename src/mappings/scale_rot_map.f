      subroutine setup_mappedgrid(rot_angle,scale)
      implicit none
      double precision rot_angle(2), scale

c     This calls maptype function that was compiled in at run time
c     It also call set_map_defaults()
c      call set_maptype()

      call set_scale(scale)
      call set_rotation(rot_angle)

      end

      subroutine set_map_defaults()
      implicit none
      double precision rot_angle(2), scale

      rot_angle(1) = 0.d0
      rot_angle(2) = 0.d0
      call set_rotation(rot_angle)

      scale = 1.d0
      call set_scale(scale)

      end


      subroutine rotate_map(xp,yp,zp)
      implicit none

      double precision xp,yp,zp
      double precision v(3), vrot(3)
      integer i,k
      double precision rrot(3,3)

      v(1) = xp
      v(2) = yp
      v(3) = zp

      call get_rotation(rrot)

      do i = 1,3
         vrot(i) = v(i)
      enddo

c     # Rotate mapping
      do i = 1,3
         vrot(i) = 0
         do k = 1,3
            vrot(i) = vrot(i) + rrot(i,k)*v(k)
         enddo
      enddo

      xp = vrot(1)
      yp = vrot(2)
      zp = vrot(3)

      end

      subroutine scale_map(xp,yp,zp)
      implicit none

      double precision xp,yp,zp
      double precision s, get_scale

      s = get_scale()

      xp = s*xp
      yp = s*yp
      zp = s*zp

      end

      subroutine set_rotation(rot_angle)
      implicit none

      double precision rot_angle(2), scale

      double precision rrot(3,3), r1(3,3), r2(3,3)

      double precision th, phi
      integer i,j,k

      logical isflat

      double precision rrot_com(3,3)
      common /comrot/ rrot_com

c     # Rotates map so as to not bias the solution.
      do i = 1,3
         do j = 1,3
            r1(i,j) = 0
         enddo
         r1(i,i) = 1
      enddo

      th = rot_angle(1)
      r1(1,1) = cos(th)
      r1(1,2) = -sin(th)
      r1(2,1) = -r1(1,2)
      r1(2,2) = r1(1,1)

      if (isflat()) then
c        # rotate in one direction only
         do i = 1,3
            do j = 1,3
               rrot(i,j) = r1(i,j)
            enddo
         enddo
      else
c        # rotate in second direction
         do i = 1,3
            do j = 1,3
               r2(i,j) = 0
            enddo
            r2(i,i) = 1
         enddo
         phi = rot_angle(2)
         r2(1,1) = cos(phi)
         r2(1,3) = -sin(phi)
         r2(3,1) = -r2(1,3)
         r2(3,3) = r2(1,1)
         do i = 1,3
            do j = 1,3
               rrot(i,j) = 0
               do k = 1,3
                  rrot(i,j) = rrot(i,j) + r2(i,k)*r1(k,j)
               enddo
            enddo
         enddo
      endif

c      call set_rotation(rrot)

      do j = 1,3
         do i = 1,3
            rrot_com(i,j) = rrot(i,j)
         enddo
      enddo

c      Writing a file like this is problematic in parallel.
c      Only one rank should write this file.
c
c      open(10,file='rrot.dat')
c      do i = 1,3
c         write(10,100) (rrot(i,j),j=1,3)
c      enddo
c      close(10)
c  100 format(3F24.16)

      end

      subroutine get_rotation(rrot)
      implicit none

      double precision rrot(3,3), rrot_com(3,3)
      common /comrot/ rrot_com

      integer i,j

      do j = 1,3
         do i = 1,3
            rrot(i,j) = rrot_com(i,j)
         enddo
      enddo

      end

      subroutine set_scale(scale)
      implicit none
      double precision scale, scale_com

      common /comscale/ scale_com
      scale_com = scale

      end


      double precision function get_scale()
      implicit none
      double precision scale_com
      common /comscale/ scale_com

      get_scale = scale_com

      end
