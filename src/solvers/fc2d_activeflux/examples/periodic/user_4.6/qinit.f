c     =====================================================
      subroutine activeflux_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &     xlower,ylower, dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

      double precision a(3,3)

      double precision pi,  pi2
      common /compi/ pi, pi2


      f(x,y)=dsin(8d0*datan(1d0)*x)*dsin(8d0*datan(1d0)*y)


      a(1,1)=1d0
      a(1,2)=4d0
      a(1,3)=1d0
      a(2,1)=4d0
      a(2,2)=16d0
      a(2,3)=4d0
      a(3,1)=1d0
      a(3,2)=4d0
      a(3,3)=1d0

      do 20 j=1-mbc,my+mbc
         do 10 i=1-mbc,mx+mbc
            q(i,j,1)=0
   10       continue
   20    continue

      do 22 j=1-mbc,my+mbc
         yj = ylower + (j-1d0)*dy
         do 12 i=1-mbc,mx+mbc
            xi = xlower + (i-1d0)*dx
            do 33 k=0,2
               do 44 l=0,2                  
                  q(i,j,1)=q(i,j,1)+(1d0/36d0)*a(k+1,l+1)
     &                     * q0(dx*k/2d0+xi,dy*l/2d0+yj)                  
   44             continue
   33          continue
   12       continue
   22    continue

      return
      end

      double precision function q0(x,y)
      implicit none

      double precision x,y

      integer example
      common /example_comm/ example  

      double precision r, r0, x0, y0, Hsmooth

 
      if (example .eq. 0) then
          q0 = dsin(8d0*datan(1d0)*x)*dsin(8d0*datan(1d0)*y)
      elseif (example .eq. 1) then
          r = sqrt((x-0.5)**2 + (y-0.5)**2)
          if (r .le. 0.2) then
              q0 = 1
          else
              q0 = 0
          endif
      elseif (example .eq. 2) then
          x0 = 0.5
          y0 = 0.5
          r0 = 0.15d0
          r = sqrt((x - x0)**2 + (y-y0)**2)
c          q0 = Hsmooth(r + r0) - Hsmooth(r - r0)
          q0 = 1 - Hsmooth(r - r0)
      endif


      end

      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end



