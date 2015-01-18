      subroutine mountain_setup()
      implicit none

      integer i
      double precision mount(469,2)

      common /mountcom/ mount

      open(10,file='mountain.dat')

      do i = 1,469
         read(10,*) mount(i,1),mount(i,2)
      enddo
      close(10)


      end

      double precision function mountain_height(xp)
      implicit none
      double precision xp, zi, zip1, xi, xip1,zm
      integer ii

      double precision mount(469,2)
      common /mountcom/ mount

      do ii = 1,468
         zi = mount(ii,2)
         zip1 = mount(ii+1,2)
         xi = mount(ii,1)
         xip1 = mount(ii+1,1)
         if (xi .le. xp .and. xp .le. xip1) then
            zm = zi + ((zip1-zi)/(xip1-xi))*(xp-xi)
            exit
         endif
      enddo

      mountain_height = zm

      end

      subroutine mapc2m_mountain(blockno,xc,yc,xp,yp,zp,scale)
      implicit none

      double precision xc,yc,xp,yp,zp, scale(2)
      integer blockno
      double precision zm, mountain_height

c     # This assumes that we have [0,1]x[0,1] coordinates
c     # coming in.

      xp = scale(1)*xc
      yp = scale(2)*yc

      zm = mountain_height(xp)

      yp = zm + (yp/scale(2))*(scale(2) - zm)

      zp = 0


      end
