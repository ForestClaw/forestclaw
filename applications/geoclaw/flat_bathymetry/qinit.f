c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise
c
       use geoclaw_module, only: sea_level, coordinate_system
       use geoclaw_module, only: grav


       implicit none
       integer i,j,meqn,mbc,mx,my,maux
       double precision xlower,ylower,dx,dy,xi,yj
       double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
       double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
c
       q = 0.d0
       do 20 i=1-mbc,mx+mbc
c         # (xi, yj) is the locaiton of the physical domain
          xi = xlower + (i-0.5d0)*dx
          do 20 j=1-mbc,my+mbc
             yj = ylower + (j-0.5d0)*dy
             if (xi.lt.9.0 .and. xi.gt. -9.0) then
c
c               # Right-going wave (square wave)
c               # q0(1,:,:) = I_([-9,9]x[-100,100])
c               # q0(2,:,:) = (u0+sqrt(g*h0)) * I_([-95,-85]x[-100,100])
c
c                q(1,i,j) = 1
c                q(2,i,j) = sqrt(grav*100) * q(1,i,j)
c
c               # Both directions (Gaussian)               
c
                q(1,i,j) = sea_level - aux(1,i,j) + 
     &                     exp(-xi**2/18.0)
             else
c
c               # Right-going wave (square wave)
c
c                q(1,i,j) = sea_level
c
c               # Both directions (Gaussian)               
c
                q(1,i,j) = sea_level - aux(1,i,j)               
             endif
  20         continue
       return
       end
