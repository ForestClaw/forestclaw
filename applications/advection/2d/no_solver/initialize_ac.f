       subroutine initialize_ac(mx,my,meqn,mbc,xlower,ylower,dx,dy,
     &      q)
       implicit none

       integer meqn,mbc,mx,my, maux
       double precision xlower, ylower,dx,dy

       double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
       double precision rp, xp, yp

       integer i,j,m
       double precision phi, x0, y0, th, pi

       pi = 4.d0*atan(1.d0)

       do j = 1-mbc,my+mbc
          do i = 1-mbc,my+mbc
             xp = xlower + (i-0.5)*dx
             yp = ylower + (j-0.5)*dy
             q(i,j,1) = -1
c            # map [0,1]x[0,1] into [-1,1]x[-1,1]
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
