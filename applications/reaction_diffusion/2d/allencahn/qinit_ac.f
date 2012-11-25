       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
       implicit none

       integer maxmx, maxmy, meqn,mbc,mx,my, maux
       double precision xlower, ylower,dx,dy

       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
       double precision xc,yc,xp,yp,zp,rp

       integer i,j,m
       double precision phi, x0, y0, th, pi

       common /compi/pi

       do j = 1-mbc,my+mbc
          do i = 1-mbc,my+mbc
             xc = xlower + (i-0.5)*dx
             yc = ylower + (j-0.5)*dy
             q(i,j,1) = -1
c            # map [0,1]x[0,1] into [-1,1]x[-1,1]
             call mapc2m(xc,yc,xp,yp,zp)
             rp = sqrt((xp)**2 + (yp)**2)
             if (rp .lt. 0.4d0) then
                q(i,j,1) = 1.d0
             endif

             do m = 1,4
c               # Larger circles
                th = (m-1)*pi/2.d0
                x0 = 0.6*cos(th)
                y0 = 0.6*sin(th)
                rp = sqrt((xp - x0)**2 + (yp-y0)**2)
                if (rp .lt. 0.2) then
                   q(i,j,1) = 1.d0
                endif

c               # Smaller circles
                th = pi/4.d0 + (m-1)*pi/2.d0
                x0 = 0.55*cos(th)
                y0 = 0.55*sin(th)
                rp = sqrt((xp - x0)**2 + (yp-y0)**2)
                if (rp .lt. 0.15) then
                   q(i,j,1) = 1.d0
                endif
             enddo
          enddo
       enddo

       return
       end
