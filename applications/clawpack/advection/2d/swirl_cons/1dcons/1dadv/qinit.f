c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, maux)

      pi = 4.d0*datan(1.d0)
c
      do i=1-mbc,mx+mbc
          q(i,1) = 1.0
      enddo
c
      return
      end
