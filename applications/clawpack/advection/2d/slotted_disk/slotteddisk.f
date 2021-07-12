      double precision function slotted_disk_sum(x,y,z)
      implicit none

      double precision x,y,z

      double precision wcloc(3,2)
      common /location_parms/ wcloc

      double precision r, hmax, b, c
      common /slotteddisk_parms/ r, hmax, b, c

      integer m, k
      double precision w(3), q, qv, slotted_disk

      q = 0
      do m = 1,2
         do k = 1,3
            w(k) = wcloc(k,m)
         end do
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
      double precision get_gc_distance, ri

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

