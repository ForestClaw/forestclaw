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

      double precision xp, yp, zp, q

      integer m, get_n_locations, n
      double precision w(3), qv, a, hmax
      double precision gaussian_flat, gaussian_sphere
      logical isflat

      n = get_n_locations()

c     # Get time dependent (td) parameters
      call get_td_gaussian_parms(a,hmax)

      q = 0
      do m = 1,n
         call get_td_gaussian_locations(m,w)
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

      double precision x,y,z,w(3),q, a, r2, hmax

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


      integer function get_n_locations()
      implicit none

      integer n_com
      common /comexp2/ n_com

      get_n_locations = n_com


      end

      subroutine set_n_locations(n)
      implicit none

      integer n

      integer n_com
      common /comexp2/ n_com

      n_com = n

      end


      subroutine set_initial_gaussian_locations(w,n)
      implicit none

      integer i, n, m
      double precision w(3,n)

      double precision wc_init_com(3,100), a_init_com, hmax_init_com
      common /comexpinit1/ wc_init_com, a_init_com, hmax_init_com

c     # Locations of initial exponentials
      if (n .gt. 100) then
         write(6,*) 'set_initial_gaussian_locations : n > 100'
         stop
      endif

      do i = 1,n
         do m = 1,3
            wc_init_com(m,i) = w(m,i)
         enddo
      enddo

      call set_n_locations(n)
      call set_td_gaussian_locations(w,n)

      end

      subroutine get_initial_gaussian_locations(i,w)
      implicit none

      integer i,m
      double precision w(3)

      double precision wc_init_com(3,100), a_init_com, hmax_init_com

      integer n, get_n_locations
      common /comexpinit1/ wc_init_com, a_init_com, hmax_init_com

      n  = get_n_locations()

      if (i .gt. n) then
         write(6,*) 'get_initial_gaussian_locations : i .gt. n'
         stop
      endif

      do m = 1,3
         w(m) = wc_init_com(m,i)
      enddo

      end


c     # ------------------------------------------
c     # Set time dependent Gaussian parms
c     # ------------------------------------------
      subroutine set_td_gaussian_locations(w,n)
      implicit none

      integer i, n, m
      double precision w(3,n)

      double precision wc_td_com(3,100), a_td_com
      double precision hmax_td_com

      common /comexptd1/ wc_td_com, a_td_com, hmax_td_com

      if (n .gt. 100) then
         write(6,*) 'set_td_gaussian_locations : n > 100'
         stop
      endif

c     # Locations of initial exponentials
      do i = 1,n
         do m = 1,3
            wc_td_com(m,i) = w(m,i)
         enddo
      enddo
      end


      subroutine get_td_gaussian_locations(i,w)
      implicit none

      integer i
      double precision w(3)

      integer m, n, get_n_locations

      double precision wc_td_com(3,100), a_td_com
      double precision hmax_td_com
      integer n_com

      common /comexptd1/ wc_td_com, a_td_com, hmax_td_com

      n = get_n_locations()
      if (i .gt. n) then
         write(6,*) 'get_td_gaussian_locations : i .gt. n'
         stop
      endif

      do m = 1,3
         w(m) = wc_td_com(m,i)
      enddo

      end

c     # ----------------------------------------------------
c     # Set a,hmax parameters for Gaussians
c     # ----------------------------------------------------

      subroutine set_initial_gaussian_parms(a,hmax)
      implicit none

      double precision a,hmax

      double precision wc_init_com(3,100), a_init_com, hmax_init_com

      common /comexpinit1/ wc_init_com, a_init_com, hmax_init_com

      a_init_com = a
      hmax_init_com = hmax

      call set_td_gaussian_parms(a,hmax)

      end


      subroutine get_initial_gaussian_parms(a,hmax)
      implicit none

      double precision a,hmax

      double precision wc_init_com(3,100), a_init_com, hmax_init_com

      common /comexpinit1/ wc_init_com, a_init_com, hmax_init_com

      a = a_init_com
      hmax = hmax_init_com

      end

      subroutine set_td_gaussian_parms(a,hmax)
      implicit none

      double precision a,hmax

      double precision wc_td_com(3,100), a_td_com, hmax_td_com

      common /comexpinit1/ wc_td_com, a_td_com, hmax_td_com

      a_td_com = a
      hmax_td_com = hmax

      end


      subroutine get_td_gaussian_parms(a,hmax)
      implicit none

      double precision a,hmax

      double precision wc_td_com(3,100), a_td_com, hmax_td_com

      common /comexpinit1/ wc_td_com, a_td_com, hmax_td_com

      a = a_td_com
      hmax = hmax_td_com

      end


c     # ------------------------------------------------
c     # Great circle distance on sphere
c     # ------------------------------------------------
      double precision function get_gc_distance(x,y,z,w)
      implicit none

      double precision x,y,z,w(3)
      double precision p(3), pn, wn, ca, th
      double precision rsphere, get_scale
      integer m

c      rsphere = get_scale()
      rsphere = 1.0

      p(1) = x
      p(2) = y
      p(3) = z
      pn = sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3))
      if (abs(pn - rsphere) .ge. 1e-3) then
         write(6,*) 'get_gc_distance : Point is not on the sphere; '
         write(6,*) x, y, z, abs(pn-rsphere)
         stop
      endif

      wn = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
      do m = 1,3
         w(m) = w(m)/wn
      enddo
      ca = (w(1)*x + w(2)*y + w(3)*z)/pn
      if (abs(ca) .gt. 1) then
         write(6,*) 'get_gc_distance : abs(ca) > 1; ', ca
         write(6,*) w(1), w(2), w(3), x,y,z
         stop
      endif
      th = acos(ca)
      get_gc_distance = th*rsphere

      end
