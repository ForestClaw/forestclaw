

c
c
c =========================================================
      subroutine src1(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
c =========================================================
      implicit none

      integer maxmx,meqn,mbc,mx,maux
      double precision t,dt,xlower,dx
      double precision q(1-mbc:maxmx+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, maux)

      double precision xcell,ux,F1,F2,F3,F4
      integer i
c
c     # 4th-order Runge-Kutta method to integrate q_t = -u_x*q
c
      do i = 1,mx
         xcell = xlower + (i-0.5d0)*dx
         ux = aux(i,2)
         F1 = -ux*q(i,1)
         F2 = -ux*(q(i,1) + 0.5d0*dt*F1)
         F3 = -ux*(q(i,1) + 0.5d0*dt*F2)
         F4 = -ux*(q(i,1) + dt*F3)
         q(i,1) = q(i,1) + dt*(F1 + 4*(F2+F3)/2.d0 + F4)/6.d0
      enddo
c
      return
      end
