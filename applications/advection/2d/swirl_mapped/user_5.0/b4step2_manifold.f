      subroutine b4step2_manifold(mx,my,mbc,dx,dy,blockno,
     &      xd,yd,zd,t, dt,maux,aux)
      implicit none

      integer mbc, mx, my, maux
      integer blockno
      double precision xlower, ylower, dx, dy, dt, t
      double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i, j
      double precision tperiod, pi2, vt, xll,yll, psi

      common /comvt/ tperiod,pi2

      include 'metric_terms.i'
c
      if (tperiod .eq. 0.d0) then
c        # special case --- indication that velocities specified in
c        # setaux should be used for all time.
         return
      endif

      call compute_velocity_psi(mx,my,mbc,
     &      dx,dy,blockno,t,xd,yd,zd,aux,maux)

      vt = cos(pi2*(t+dt/2.d0)/tperiod)


      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # multiply by time-factor:
            aux(1,i,j) = vt * aux(1,i,j)
            aux(2,i,j) = vt * aux(2,i,j)
         enddo
      enddo

      return
      end
