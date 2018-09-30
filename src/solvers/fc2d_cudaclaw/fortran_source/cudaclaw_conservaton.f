c
c
c     =====================================================
      subroutine cudaclaw_flux_add(mx,my,mbc,meqn,dx,dy,
     &      dt, qnew,flux,iface, buffer)
c     =====================================================

      implicit none

      integer mx,my,mbc,meqn, iface
      double precision dx,dy, dt
      double precision qnew(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision flux(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision buffer(max(mx,my),meqn)

      integer i,j,m

      do m = 1,meqn
         if (iface == 0) then
            do j = 1,my
               buffer(j,m) = buffer(j,m) + flux(1,j,m)
            enddo
         endif

         if (iface == 1) then
            do j = 1,my
               buffer(j,m) = buffer(j,m) + flux(mx+1,j,m)
            enddo
         endif

         if (iface == 2) then
            do i = 1,mx
               buffer(i,m) = buffer(i,m) + flux(i,1,m)
            enddo
         endif

         if (iface == 3) then
            do i = 1,mx
               buffer(i,m) = buffer(i,m) + flux(i,my+1,m)
            enddo
         endif

      enddo

      end
