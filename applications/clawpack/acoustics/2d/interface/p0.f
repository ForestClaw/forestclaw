c
c
c
c     =====================================================
      double precision function p0(r)
c     =====================================================
c
      implicit real*8(a-h,o-z)
c
c     # Initial variation in pressure with respect to r
c
c
      width = 0.1d0
      pi = 4.d0*datan(1.d0)    
      if (r .lt. width) then
              p0 = 1.d0 + 0.5d0 * (dcos(pi*r/width) - 1.d0)
         else
              p0 = 0.d0
         endif
      return
      end
