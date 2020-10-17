       subroutine allencahn_init(meqn,mbc,mx,my,
     &                   xlower,ylower,dx,dy,q)
       implicit none

       integer meqn,mbc,mx,my
       double precision xlower, ylower,dx,dy

       double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
       double precision xc,yc,xp,yp,zp,rp

       integer i,j,m
       double precision phi, x0, y0, th

       double precision pi
       common /compi/pi

       do j = 1-mbc,my+mbc
          do i = 1-mbc,mx+mbc
             xc = xlower + (i-0.5)*dx
             yc = ylower + (j-0.5)*dy
             q(i,j,1) = -1
c             call mapc2m(xc,yc,xp,yp,zp)
             xp = -1 + 2*xc
             yp = -1 + 2*yc
             zp = 0
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
                x0 = 0.55d0*cos(th)
                y0 = 0.55d0*sin(th)
                rp = sqrt((xp - x0)**2 + (yp-y0)**2)
                if (rp .lt. 0.15) then
                   q(i,j,1) = 1.d0
                endif
             enddo
          enddo
       enddo

       return
       end
