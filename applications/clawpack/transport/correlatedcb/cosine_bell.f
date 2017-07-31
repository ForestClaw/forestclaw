      double precision function cosine_bell_sum(x,y,z)
      implicit none

      double precision x,y,z

      double precision hmax, h1, b, c, r, w(3)
      double precision q, qv

      integer m,n, get_n_locations
      double precision cosbell

      n = get_n_locations()

      call get_td_cosbell_parms(r,hmax,b,c)

      q = 0
      do m = 1,n
         call get_td_gaussian_locations(m,w)
         qv = cosbell(x,y,z,w,r,hmax)
         q  = q  + qv
      enddo
      cosine_bell_sum = b + c*q

      end


      double precision function cosbell(x,y,z,w,r,hmax)
      implicit none

      double precision x,y,z,w(3),r,hmax
      double precision get_gc_distance, h1, ri
      double precision pi

      common /compi/ pi


      ri = get_gc_distance(x,y,z,w)

      if (ri .le. r) then
         h1 = hmax*(1 + cos(pi*ri/r))
      else
         h1 = 0.d0
      endif
      cosbell = h1

      end


c     # ----------------------------------------------------
c     # Set a,hmax parameters for Cosine bell
c     # ----------------------------------------------------

      subroutine set_initial_cosbell_parms(r,hmax,b,c)
      implicit none

      double precision r,hmax, b,c

      double precision wc_init_com(3,100), r_init_com, hmax_init_com
      double precision b_init_com, c_init_com

      common /comcbinit1/ wc_init_com, r_init_com, hmax_init_com,
     &      b_init_com, c_init_com

      r_init_com = r
      hmax_init_com = hmax
      b_init_com = b
      c_init_com = c

      call set_td_cosbell_parms(r,hmax,b,c)

      end


      subroutine get_initial_cosbell_parms(r,hmax,b,c)
      implicit none

      double precision r,hmax, b,c

      double precision wc_init_com(3,100), r_init_com,
     &      hmax_init_com, b_init_com, c_init_com

      common /comcbinit1/ wc_init_com, r_init_com, hmax_init_com,
     &      b_init_com, c_init_com

      r = r_init_com
      hmax = hmax_init_com
      b = b_init_com
      c = c_init_com

      end

      subroutine set_td_cosbell_parms(r,hmax,b,c)
      implicit none

      double precision r,hmax,b,c

      double precision wc_td_com(3,100), r_td_com, hmax_td_com,
     &      b_td_com, c_td_com

      common /comcbinit1/ wc_td_com, r_td_com, hmax_td_com,
     &      b_td_com, c_td_com

      r_td_com = r
      hmax_td_com = hmax
      b_td_com = b
      c_td_com = c

      end


      subroutine get_td_cosbell_parms(r,hmax,b,c)
      implicit none

      double precision r,hmax, b, c

      double precision wc_td_com(3,100), r_td_com, hmax_td_com
      double precision b_td_com,c_td_com

      common /comcbinit1/ wc_td_com, r_td_com, hmax_td_com,
     &      b_td_com, c_td_com

      r = r_td_com
      hmax = hmax_td_com
      b = b_td_com
      c = c_td_com

      end
