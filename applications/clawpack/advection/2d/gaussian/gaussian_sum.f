c     # -----------------------------------------------------
c     # GAUSSIAN_SUM
c     #
c     # Computes the value q at a point (x,y,z), where q
c     # is defined by the sum
c     #
c     #             q = sum_{i=1,n} g(x,y,z;xi,yi,zi,a,hmax)
c     #
c     #          g(x,y,z;xi,yi,zi,a,hmax) = hmax*exp(-a*r2)
c     #
c     #          r2 = (x-xi)**2 + (y-yi)**2 + (z-zi)**2
c     #
c     # Note : 'td' means "time dependent"
c     # -----------------------------------------------------

      double precision function gaussian_sum(xp,yp,zp)
      implicit none

      double precision xp, yp, zp

      double precision wc(3,2)
      common /location_parms/ wc

      double precision a, hmax
      common /gaussian_parms/ a, hmax

      integer m, k
      double precision w(3), qv, q
      double precision gaussian_flat, gaussian_sphere
      logical isflat

      q = 0
      do m = 1,2
         do k = 1,3
            w(k) = wc(k,m)
         end do
         if (isflat()) then
            qv = gaussian_flat(xp,yp,zp,w,a,hmax)
         else
            qv = gaussian_sphere(xp,yp,zp,w,a,hmax)
         endif
         q  = q  + qv
      enddo

      gaussian_sum = q

      end

      double precision function gaussian_flat(x,y,z,w,a,hmax)
      implicit none

      double precision x,y,z,w(3),a, r2, hmax

      r2 = (x-w(1))**2 + (y - w(2))**2 + (z - w(3))**2
      gaussian_flat = hmax*exp(-a*r2)

      end


      double precision function gaussian_sphere(x,y,z,w,a,hmax)
      implicit none

      double precision x,y,z,w(3),a,q, hmax

      double precision get_gc_distance, d

      d = get_gc_distance(x,y,z,w)
      q = hmax*exp(-a*d**2)

      gaussian_sphere = q

      end
