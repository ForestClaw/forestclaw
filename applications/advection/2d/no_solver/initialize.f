      subroutine initialize(mx,my,meqn,mbc,xlower,ylower,dx,dy,q)
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my
      double precision xlower, ylower, dx, dy
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)

      integer i, j, mq, ichoice
      double precision xlow, ylow,wl

      common /com_init/ ichoice

      if (meqn .gt. 2) then
         write(6,*) 'initialize (initialize.f) : meqn > 2'
         stop
      endif

      ichoice = 1


      do i = 1-mbc,mx+mbc
         xlow = xlower + (i-1)*dx
         do j = 1-mbc,my+mbc
            ylow = ylower + (j-1)*dy

            call cellave2(xlow,ylow,dx,dy,wl)

            q(i,j,1) = wl
            if (meqn > 1) then
               q(i,j,2) = 1-q(i,j,1)
            endif
         enddo
      enddo

      return
      end
