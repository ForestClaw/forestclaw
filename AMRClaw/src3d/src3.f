      subroutine src3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &      xlower, ylower, zlower, dx, dy, dz, q,maux,aux,t,dt)
      implicit none

      integer maxmx, maxmy, maxmz, meqn, mbc, mx, my, mz,maux
      double precision xlower, ylower, zlower, dx, dy, dz, t,dt

      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux)

c
c      # dummy subroutine for use when equation has no source term.
c      # If method(5)=0 then this routine is never called, but its
c      # existence may be required by some compilers.
c
      return
      end
