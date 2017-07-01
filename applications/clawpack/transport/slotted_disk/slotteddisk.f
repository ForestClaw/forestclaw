      double precision function slotted_disk_sum(x,y,z)
      implicit none

      double precision x,y,z

      double precision w(3), q, qv, b, c
      double precision r, hmax
      integer m, n, get_n_locations
      double precision slotted_disk

      n = get_n_locations()

      call get_td_sdisk_parms(r,hmax,b,c)

      q = 0
      do m = 1,n
         call get_td_gaussian_locations(m,w)
         qv = slotted_disk(x,y,z,w,r,hmax,m)
         q = q + qv
      enddo

      slotted_disk_sum = b + c*q
      end

      double precision function slotted_disk(x,y,z,w,r,hmax,
     &      which_disk)
      implicit none

      double precision x,y,z,w(3),r,hmax
      integer which_disk

      double precision th, lambda, thi, li, q
      double precision get_gc_distance, ri, pi

c     r = 0.5d0

      ri = get_gc_distance(x,y,z,w)

      call map2polar(w(1),w(2),w(3),lambda,th)
      call map2polar(x,y,z,li,thi)

c     # Default value
      q = 0

c     # Now check to if we are in a disk
      if (ri .lt. r) then
         if (abs(lambda - li) .gt. r/6.d0) then
            q = 1.d0
         else
            if (which_disk .eq. 1) then
               if (th - thi .gt. 5.d0*r/12) then
                  q = 1.d0
               endif
            else
               if (th - thi .lt. -5.d0*r/12) then
                  q = 1.d0
               endif
            endif
         endif
      endif

      slotted_disk = q

      end


c     # ----------------------------------------------------
c     # Set a,hmax parameters for slotted disk
c     # ----------------------------------------------------

      subroutine set_initial_sdisk_parms(r,hmax,b,c)
      implicit none

      double precision r,hmax, b,c

      double precision wc_init_com(3,100), r_init_com, hmax_init_com
      double precision b_init_com, c_init_com

      common /comsdinit1/ wc_init_com, r_init_com, hmax_init_com,
     &      b_init_com, c_init_com

      r_init_com = r
      hmax_init_com = hmax
      b_init_com = b
      c_init_com = c

      call set_td_sdisk_parms(r,hmax,b,c)

      end


      subroutine get_initial_sdisk_parms(r,hmax,b,c)
      implicit none

      double precision r,hmax, b,c

      double precision wc_init_com(3,100), r_init_com,
     &      hmax_init_com, b_init_com, c_init_com

      common /comsdinit1/ wc_init_com, r_init_com, hmax_init_com,
     &      b_init_com, c_init_com

      r = r_init_com
      hmax = hmax_init_com
      b = b_init_com
      c = c_init_com

      end

      subroutine set_td_sdisk_parms(r,hmax,b,c)
      implicit none

      double precision r,hmax,b,c

      double precision wc_td_com(3,100), r_td_com, hmax_td_com,
     &      b_td_com, c_td_com

      common /comsdinit1/ wc_td_com, r_td_com, hmax_td_com,
     &      b_td_com, c_td_com

      r_td_com = r
      hmax_td_com = hmax
      b_td_com = b
      c_td_com = c

      end


      subroutine get_td_sdisk_parms(r,hmax,b,c)
      implicit none

      double precision r,hmax, b, c

      double precision wc_td_com(3,100), r_td_com, hmax_td_com
      double precision b_td_com,c_td_com

      common /comsdinit1/ wc_td_com, r_td_com, hmax_td_com,
     &      b_td_com, c_td_com

      r = r_td_com
      hmax = hmax_td_com
      b = b_td_com
      c = c_td_com

      end
