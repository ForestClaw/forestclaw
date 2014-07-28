      subroutine setup_map()
      implicit none

      integer smallend(2),bigend(2)
      integer mgrid, ngrid, ntotal
      double precision r1
      integer n, m

c     # Number of boxes
      call create_boxes()

      open(10,file='boxes.dat')
      read(10,*) ntotal
      read(10,*) mgrid
c      read(10,*) ngrid
c      read(10,*) r1
      do n = 1,ntotal
         read(10,*) smallend(1),smallend(2), ngrid, r1
         do m = 1,2
            bigend(m) = smallend(m) + ngrid
         enddo
         call add_box(Mgrid,smallend,bigend,r1,ngrid,n)
      enddo
      close(10)
      end

      subroutine create_boxes()
      implicit none

      integer n_box_com

      common /combox0/ n_box_com
      n_box_com = 0

      end

      subroutine add_box(mgrid,smallend,bigend,r1,ngrid,i)
      implicit none

      integer Mgrid, smallend(2), bigend(2), ngrid,i
      double precision boxes(5), r1
      double precision xlow,ylow,xhi,yhi,h

      double precision boxes_com(6,100)
      integer n_box_com

      common /combox0/ n_box_com
      common /combox1/ boxes_com

      n_box_com = n_box_com + 1

      h = 2.d0/mgrid
      xlow = -1 + (smallend(1)-1)*h
      xhi = -1 + (bigend(1)-1)*h
      ylow = -1 + (smallend(2)-1)*h
      yhi = -1 + (bigend(2)-1)*h

      boxes_com(1,i) = xlow
      boxes_com(2,i) = ylow
      boxes_com(3,i) = xhi
      boxes_com(4,i) = yhi
      boxes_com(5,i) = r1
      boxes_com(6,i) = ngrid

      end


      subroutine check_boxes(dx,dy)
      implicit none

      double precision dx,dy

      integer n, m
      logical check(4)

      double precision xlow,ylow,xhi,yhi,r1,ts
      double precision ndx,ndy, rx,ry
      double precision ndx1,ndy1, rx1, ry1

      double precision boxes_com(6,100)
      integer n_box_com

      common /combox0/ n_box_com
      common /combox1/ boxes_com

c     # Only print out information for one box, since we assume
c     # all boxes are the same size now.
      do n = 1,n_box_com
         xlow = boxes_com(1,n)
         ylow = boxes_com(2,n)
         xhi  = boxes_com(3,n)
         yhi  = boxes_com(4,n)
         r1   = boxes_com(5,n)
c         ngrid = boxes_com(6,n)

         ndx = (xhi - xlow)/dx
         ndy = (yhi - ylow)/dy
         ndx1 = nint(ndx*ts)/ts
         ndy1 = nint(ndy*ts)/ts

         ts = 1.d5
         rx = r1*ndx
         ry = r1*ndy
         rx1 = nint(rx*ts)/ts
         ry1 = nint(ry*ts)/ts

         check(1) = abs(ndx - ndx1) .lt. 1d-8
         check(2) = abs(ndy - ndy1) .lt. 1d-8
         check(3) = abs(rx - rx1) .lt. 1d-8
         check(4) = abs(ry - ry1) .lt. 1d-8

         write(6,*) ' '
         write(6,95) '------------'
         write(6,90) 'Box ', n
         write(6,95) '------------'
   90    format(A,I2)
   95    format(A)
         if (check(1)) then
            write(6,100) 'Square : x direction',nint(ndx)
         else
            write(6,110) 'Square : x direction',ndx
         endif

         if (check(2)) then
            write(6,100) 'Square : y direction',nint(ndy)
         else
            write(6,110) 'Square : y direction',ndy
         endif

         if (check(3)) then
            write(6,100) 'Circle : x direction',nint(rx)
         else
            write(6,110) 'Circle : x direction',rx
         endif

         if (check(4)) then
            write(6,100) 'Circle : y direction',nint(ry)
         else
            write(6,110) 'Circle : y direction',ry
         endif

         do m = 1,4
            if (.not. check(m)) then
               write(6,*) '*** WARNING : Boxes or circles do not ',
     &               ' align with grid.'
               stop
            endif
         enddo
      enddo
      write(6,*) ' '

  100 format(A20,I10)
  110 format(A20,F10.5)
      end

      logical function is_in_circle(xp,yp,dr)
      implicit none

      double precision xp,yp, xc1, yc1, zp

      double precision boxes_com(6,100)
      integer n_box_com

      common /combox0/ n_box_com
      common /combox1/ boxes_com

      integer i
      double precision r, xc, yc, r2, dr

      is_in_circle = .false.
      do i = 1,n_box_com
         xc = (boxes_com(1,i) + boxes_com(3,i))/2
         yc = (boxes_com(2,i) + boxes_com(4,i))/2

         r = boxes_com(5,i)*(boxes_com(3,i) - boxes_com(1,i))/2

         r2 = sqrt((xp-xc)**2 + (yp - yc)**2)
         if (abs(r2 - r)/r .le. 0.1d0) then
            is_in_circle = .true.
            return
         endif
      enddo

      end
