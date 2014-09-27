c     # ----------------------------------------------
c     # Scaling routines
c     # ----------------------------------------------

c      subroutine set_scale(scale)
c      implicit none
c      double precision scale(3), scale_com(3)
c      integer m
c
c      common /comscale/ scale_com
c
c      do m = 1,3
c         scale_com(m) = scale(m)
c      enddo
c
c      end
c
c
c      subroutine  get_scale(scale)
c      implicit none
c      double precision scale_com(3), scale(3)
c      integer m
c      common /comscale/ scale_com
c
c      do m = 1,3
c         scale(m) = scale_com(m)
c      enddo
c
c
c      end
c
c      subroutine scale_map(xp,yp,zp)
c      implicit none
c
c      double precision xp,yp,zp
c      double precision s(3)
c
c      call get_scale(scale)
c
c      xp = s(1)*xp
c      yp = s(2)*yp
c      zp = s(3)*zp
c
c      end


c     # ----------------------------------------------
c     # Shift routines
c     # ----------------------------------------------
c
c      subroutine set_shift(shift)
c      implicit none
c
c      double precision shift(3), shift_com(3)
c      integer m
c
c      common /comshift/ shift_com
c
c      do m = 1,3
c         shift_com(m) = shift(m)
c      enddo
c
c
c      end
c
c      subroutine get_shift(shift)
c      implicit none
c
c      double precision shift(3), shift_com(3)
c      integer m
c
c      common /comshift/ shift_com
c
c      do m = 1,3
c         shift(m) = shift_com(m)
c      enddo
c
c
c      end
c
c      subroutine shift_map(xp,yp,zp)
c      implicit none
c
c      double precision xp,yp,zp
c      double precision shift(3)
c      integer m
c
c
c      call get_shift(shift)
c
c      xp = xp + shift(1)
c      yp = yp + shift(2)
c      zp = zp + shift(3)
c      end

c     # ----------------------------------------------
c     # Rotation routines
c     # ----------------------------------------------
      subroutine set_rotation_matrix(rot_angle,rrot)
      implicit none

      double precision rot_angle(2)
      double precision rrot(3,3)

      double precision r1(3,3), r2(3,3)

      double precision th, phi
      integer i,j,k

c      double precision rrot_com(3,3)
c      common /comrot/ rrot_com

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

c     # rotate in second direction (for full 3d map)
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

c      do j = 1,3
c         do i = 1,3
c            rrot_com(i,j) = rrot(i,j)
c         enddo
c      enddo

      end


c      subroutine get_rotation(rrot)
c      implicit none
c
c      double precision rrot(3,3), rrot_com(3,3)
c      common /comrot/ rrot_com
c
c      integer i,j
c
c      do j = 1,3
c         do i = 1,3
c            rrot(i,j) = rrot_com(i,j)
c         enddo
c      enddo
c
c      end


c      subroutine rotate_mapping(rrot,xp,yp,zp)
c      implicit none
c
c      double precision xp,yp,zp
c      double precision v(3), vrot(3)
c      integer i,k
c      double precision rrot(3,3)
c
c      v(1) = xp
c      v(2) = yp
c      v(3) = zp
c
c      call get_rotation(rrot)
c
c      do i = 1,3
c         vrot(i) = v(i)
c      enddo
c
cc     # Rotate mapping
c      do i = 1,3
c         vrot(i) = 0
c         do k = 1,3
c            vrot(i) = vrot(i) + rrot(i,k)*v(k)
c         enddo
c      enddo
c
c      xp = vrot(1)
c      yp = vrot(2)
c      zp = vrot(3)
c
c      end
