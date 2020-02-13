       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

      integer i,j
      double precision xc,yc

      integer init_choice
      common /initchoice_comm/ init_choice

      do j = 1-mbc,my+mbc
          do i = 1-mbc,mx+mbc     
              xc = xlower + (i-0.5)*dx
              yc = ylower + (j-0.5)*dy
              if (init_choice .eq. 0) then
                  if (abs(xc-0.5) .le. 0.25 
     &                    .and. abs(yc-0.5) .le. 0.25) then
                      q(1,i,j) = 1.d0
                  else
                      q(1,i,j) = 0.d0
                  endif
              elseif (init_choice .eq. 1) then
                  q(1,i,j) = 1.d0
              endif
          enddo
      enddo

      return
      end
